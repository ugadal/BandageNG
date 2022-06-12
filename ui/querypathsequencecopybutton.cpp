#include "querypathsequencecopybutton.h"

#include <QApplication>
#include <QClipboard>
#include <utility>

QueryPathSequenceCopyButton::QueryPathSequenceCopyButton(QByteArray pathSequence, const QString& pathStart) :
    m_pathSequence(std::move(pathSequence))
{
    setText(pathStart);
    connect(this, SIGNAL(clicked(bool)), this, SLOT(copySequenceToClipboard()));
}


void QueryPathSequenceCopyButton::copySequenceToClipboard()
{
    QClipboard * clipboard = QApplication::clipboard();
    clipboard->setText(m_pathSequence);
}
