#pragma once

#include "assemblygraph.h"
#include <QString>

class AssemblyGraphBuilder {
  public:    
    virtual void build(AssemblyGraph &graph) = 0;
    virtual ~AssemblyGraphBuilder() = default;

    static std::unique_ptr<AssemblyGraphBuilder> get(const QString &fullFileName);

    bool hasCustomLables() const { return hasCustomLabels_; }
    bool hasCustomColours() const { return hasCustomColours_; }
    bool hasComplexOverlaps() const { return hasComplexOverlaps_; }
    
  protected:
    AssemblyGraphBuilder(const QString &fileName)
            : fileName_(fileName) {}

    QString fileName_;
    bool hasCustomLabels_ = false;
    bool hasCustomColours_ = false;
    bool hasComplexOverlaps_ = false;
};




