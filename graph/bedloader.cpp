// Copyright 2022 Andrey Zakharov

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

#include "bedloader.hpp"

#include <sstream>

#include <fast-cpp-csv-parser/csv.h>


namespace bed {

std::vector<int64_t> parseIntArray(const std::string &intArrayString) {
    std::stringstream ss{intArrayString};
    std::string intString;
    std::vector<int64_t> res;
    while (std::getline(ss, intString, ',')) {
        if (intString.empty()) {
            break;
        }
        res.push_back(std::stoll(intString));
    }
    return res;
}

/* TODO: 1) add assertions
 *       2) support custom fields ?
 *       3) make better error handling
 */
std::vector<Line> load(const std::filesystem::path &path) {
    io::CSVReader<12,
                  io::trim_chars<' '>,
                  io::no_quote_escape<'\t'>,
                  io::throw_on_overflow,
                  io::single_and_empty_line_comment<'#'>> csvReader(path);
    std::vector<Line> res;
    Line bedLine;
    char strandChar = '.';
    std::string itemRgbString;
    int64_t blockCount = 0;
    std::string blockSizesString;
    std::string blockStartsString;
    while (csvReader.read_row(bedLine.chrom,
                              bedLine.chromStart,
                              bedLine.chromEnd,
                              bedLine.name,
                              bedLine.score,
                              strandChar,
                              bedLine.thickStart,
                              bedLine.thickEnd,
                              itemRgbString,
                              blockCount,
                              blockSizesString,
                              blockStartsString)) {
        bedLine.strand = Strand{strandChar};
        if (bedLine.thickStart == -1 || bedLine.thickEnd == -1) {
            bedLine.thickStart = bedLine.chromStart;
            bedLine.thickEnd = bedLine.chromEnd;
        }
        if (itemRgbString != "0") {
            auto rgbArray = parseIntArray(itemRgbString);
            bedLine.itemRgb.r = rgbArray[0];
            bedLine.itemRgb.g = rgbArray[1];
            bedLine.itemRgb.b = rgbArray[2];
        }
        if (blockCount != 0) {
            auto blockSizes = parseIntArray(blockSizesString);
            auto blockStarts = parseIntArray(blockStartsString);
            for (int i = 0; i < blockCount; i++) {
                auto blockStart = bedLine.chromStart + blockStarts[i];
                auto blockEnd = blockStart + blockSizes[i];
                bedLine.blocks.push_back(Block{.start = blockStart, .end = blockEnd});
            }
        }
        res.emplace_back(bedLine);
    }
    return res;
}

}
