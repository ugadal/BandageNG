//
// Created by wafemand on 26/02/2022.
//

#ifndef BANDAGE_GRAPH_SEQUENCEUTILS_HPP_
#define BANDAGE_GRAPH_SEQUENCEUTILS_HPP_

#include <QByteArray>
#include "seq/sequence.hpp"

inline QByteArray sequenceToQByteArray(const Sequence &sequence) {
    auto sequence_str = sequence.str();
    return {sequence_str.c_str(), static_cast<qsizetype>(sequence_str.size())};
}

#endif //BANDAGE_GRAPH_SEQUENCEUTILS_HPP_
