// Copyright 2023 Anton Korobeynikov

// This file is part of Bandage-NG

// Bandage-NG is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Bandage-NG is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Bandage.  If not, see <http://www.gnu.org/licenses/>.

#include <QApplication>
#include <filesystem>

namespace CLI {
    class App;
}

struct QueryPathsCmd {
    std::filesystem::path m_graph;
    std::filesystem::path m_queries;
    std::string m_prefix;
    bool m_pathFasta = false;
    bool m_hitsFasta = false;
};

CLI::App *addQueryPathsSubcommand(CLI::App &app,
                                  QueryPathsCmd &cmd);
int handleQueryPathsCmd(QApplication *app, const QueryPathsCmd &cmd);
