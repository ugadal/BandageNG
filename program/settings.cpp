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


#include "settings.h"
#include "graph/nodecolorer.h"
#include <QDir>

Settings::Settings()
{
    doubleMode = false;

    nodeLengthMode = AUTO_NODE_LENGTH;
    autoNodeLengthPerMegabase = 1000.0;
    manualNodeLengthPerMegabase = FloatSetting(1000.0, 0, 1000000.0);
    meanNodeLength = 40.0;
    minTotalGraphLength = 500.0;
    graphLayoutQuality = IntSetting(2, 0, 4);
    linearLayout = false;
    minimumNodeLength = FloatSetting(5.0, 1.0, 100.0);
    edgeLength = FloatSetting(5.0, 0.1, 100.0);
    doubleModeNodeSeparation = FloatSetting(2.0, 0.0, 100.0);
    nodeSegmentLength = FloatSetting(20.0, 1.0, 1000.0);
    componentSeparation = FloatSetting(50.0, 0, 1000.0);

    averageNodeWidth = FloatSetting(5.0, 0.5, 1000.0);
    depthEffectOnWidth = FloatSetting(0.5, 0.0, 1.0);
    depthPower = FloatSetting(0.5, 0.0, 1.0);

    jumpsAsLinks = false;
    edgeWidth = FloatSetting(1.5, 0.1, 100);
    linkWidth = FloatSetting(0.5, 0.1, 100);
    outlineThickness = FloatSetting(0.0, 0.0, 100.0);
    selectionThickness = 1.0;
    arrowheadsInSingleMode = false;
    textOutlineThickness = FloatSetting(1.5, 0.0, 10.0);

    blastRainbowPartsPerQuery = 100;

    graphScope = WHOLE_GRAPH;
    nodeDistance = IntSetting(0, 0, 100);
    startingNodesExactMatch = true;
    startingNodes = "";
    blastQueryFilename = "";
    unnamedQueryDefaultName = "unnamed";

    minZoom = 0.001;
    minZoomOnGraphDraw = 0.001;
    maxZoom = 100.0;
    maxAutomaticZoom = 10.0;
    zoomFactor = 1.0005;
    textZoomScaleFactor = 0.7;

    dragStrength = 100.0;

    displayNodeCustomLabels = false;
    displayNodeNames = false;
    displayNodeLengths = false;
    displayNodeDepth = false;
    displayNodeCsvData = false;
    displayNodeCsvDataCol = 0;
    labelFont = QFont();
    textOutline = false;
    antialiasing = true;
    positionTextNodeCentre = false;

    nodeDragging = NEARBY_PIECES;

    initializeColorer(RANDOM_COLOURS);
    uniformPositiveNodeColour = QColor(178, 34, 34);
    uniformNegativeNodeColour = QColor(128, 0, 0);
    uniformNodeSpecialColour = QColor(0, 128, 0);

    randomColourPositiveOpacity = IntSetting(255, 0, 255);
    randomColourNegativeOpacity = IntSetting(255, 0, 255);
    randomColourPositiveSaturation = IntSetting(127, 0, 255);
    randomColourNegativeSaturation = IntSetting(127, 0, 255);
    randomColourPositiveLightness = IntSetting(150, 0, 255);
    randomColourNegativeLightness = IntSetting(90, 0, 255);

    edgeColour = QColor(0, 0, 0, 180);
    outlineColour = QColor(0, 0, 0);
    selectionColour = QColor(0, 0, 255);
    textColour = QColor(0, 0, 0);
    textOutlineColour = QColor(255, 255, 255);

    contiguitySearchSteps = IntSetting(15, 1, 50);
    contiguousStrandSpecificColour = QColor(0, 155, 0);
    contiguousEitherStrandColour = QColor(0, 155, 0);
    maybeContiguousColour = QColor(190, 240, 190);
    notContiguousColour = QColor(220, 220, 220);
    contiguityStartingColour = QColor(0, 255, 0);

    grayColor = QColor(220, 220, 220);

    autoDepthValue = true;
    lowDepthValue = FloatSetting(5.0, 0.0, 1000000.0);
    highDepthValue = FloatSetting(50.0, 0.0, 1000000.0);
    colorMap = ColorMap::Viridis;

    pathHighlightShadingColour = QColor(0, 0, 0, 60);
    pathHighlightOutlineColour = QColor(0, 0, 0);

    defaultCustomNodeColour = QColor(190, 190, 190);

    minAutoFindEdgeOverlap = 10;
    maxAutoFindEdgeOverlap = 200;

    maxHitsForQueryPath = IntSetting(100, 0, 1000000);
    maxQueryPathNodes = IntSetting(6, 1, 50);
    minQueryCoveredByPath = FloatSetting(0.9, 0.15, 1.0);
    minQueryCoveredByHits = FloatSetting(0.9, 0.15, 1.0, true);
    minMeanHitIdentity = FloatSetting(0.5, 0.0, 1.0, true);
    maxEValueProduct = SciNotSetting(SciNot(1.0, -10), SciNot(1.0, -999), SciNot(9.9, 1), true);
    minLengthPercentage = FloatSetting(0.95, 0.0, 10000.0, true);
    maxLengthPercentage = FloatSetting(1.05, 0.0, 10000.0, true);
    minLengthBaseDiscrepancy = IntSetting(-100, -1000000, 1000000, false);
    maxLengthBaseDiscrepancy = IntSetting(100, -1000000, 1000000, false);

    blastSearchParameters = "";

    blastAlignmentLengthFilter = IntSetting(100, 1, 1000000, false);
    blastQueryCoverageFilter = FloatSetting(50.0, 0.0, 100.0, false);
    blastIdentityFilter = FloatSetting(90.0, 0.0, 100.0, false);
    blastEValueFilter = SciNotSetting(SciNot(1.0, -10), SciNot(1.0, -999), SciNot(9.9, 1), false);
    blastBitScoreFilter = FloatSetting(1000.0, 0.0, 1000000.0, false);
    blastAnnotationGroupName = "Blast Hits";
    bedAnnotationGroupName = "BED";

    minDepthRange = FloatSetting(10.0, 0.0, 1000000.0);
    maxDepthRange = FloatSetting(100.0, 0.0, 1000000.0);

    annotationsSettings = {};
}


void Settings::initializeColorer(NodeColorScheme scheme) {
    nodeColorer = INodeColorer::create(scheme);
}
