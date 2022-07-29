//Copyright 2017 Ryan Wick

//This file is part of Bandage

//Bandage is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//Bandage is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with Bandage.  If not, see <http://www.gnu.org/licenses/>.


#include "buildblastdatabaseworker.h"
#include "blastsearch.h"

#include "graph/debruijnnode.h"
#include "graph/assemblygraph.h"

#include "program/globals.h"

#include <QProcess>
#include <QFile>
#include <QTextStream>

BuildBlastDatabaseWorker::BuildBlastDatabaseWorker(QString makeblastdbCommand)
        : m_makeblastdbCommand(makeblastdbCommand) {}

bool BuildBlastDatabaseWorker::buildBlastDatabase() {
    g_blastSearch->m_cancelBuildBlastDatabase = false;

    QFile file(g_blastSearch->m_tempDirectory.filePath("all_nodes.fasta"));
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        emit finishedBuild("Failed to open: " + file.fileName());
        return false;
    }
    QTextStream out(&file);

    for (const auto *node : g_assemblyGraph->m_deBruijnGraphNodes) {
        if (g_blastSearch->m_cancelBuildBlastDatabase) {
            emit finishedBuild("Build cancelled.");
            return false;
        }

        out << node->getFasta(true, false, false);
    }
    file.close();

    // Make sure the graph has sequences to BLAST.
    bool atLeastOneSequence = false;
    for (const auto *node : g_assemblyGraph->m_deBruijnGraphNodes) {
        if (!node->sequenceIsMissing()) {
            atLeastOneSequence = true;
            break;
        }
    }

    if (!atLeastOneSequence) {
        m_error = "Cannot build the BLAST database as this graph contains no sequences";
        emit finishedBuild(m_error);
        return false;
    }

    QStringList makeBlastdbOptions;
    makeBlastdbOptions << "-in" << (g_blastSearch->m_tempDirectory.filePath("all_nodes.fasta"))
                       << "-dbtype" << "nucl";
    
    g_blastSearch->m_makeblastdb = new QProcess();
    g_blastSearch->m_makeblastdb->start(m_makeblastdbCommand,
                                        makeBlastdbOptions);

    bool finished = g_blastSearch->m_makeblastdb->waitForFinished(-1);

    if (g_blastSearch->m_makeblastdb->exitCode() != 0 || !finished) {
        m_error = "There was a problem building the BLAST database";
        QString stdErr = g_blastSearch->m_makeblastdb->readAllStandardError();
        m_error += stdErr.isEmpty() ? "." : ":\n\n" + stdErr;
    } else if (g_blastSearch->m_cancelBuildBlastDatabase)
        m_error = "Build cancelled.";
    else
        m_error = "";

    emit finishedBuild(m_error);

    g_blastSearch->m_makeblastdb->deleteLater();
    g_blastSearch->m_makeblastdb = nullptr;
    return m_error.isEmpty();
}
