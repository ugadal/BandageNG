#include "assemblygraph.h"
#include <QMapIterator>
#include "../program/globals.h"
#include "../program/settings.h"

AssemblyGraph::AssemblyGraph()
{
    m_ogdfGraph = new ogdf::Graph();
    m_graphAttributes = new ogdf::GraphAttributes(*m_ogdfGraph, ogdf::GraphAttributes::nodeGraphics |
                                                  ogdf::GraphAttributes::edgeGraphics);
}

AssemblyGraph::~AssemblyGraph()
{

}


void AssemblyGraph::cleanUp()
{
    QMapIterator<int, DeBruijnNode*> i(m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();
        delete i.value();
    }
    m_deBruijnGraphNodes.clear();

    for (size_t i = 0; i < m_deBruijnGraphEdges.size(); ++i)
        delete m_deBruijnGraphEdges[i];
    m_deBruijnGraphEdges.clear();
}





//This function makes a double edge: in one direction for the given nodes
//and the opposite direction for their reverse complements.  It adds the
//new edges to the vector here and to the nodes themselves.
void AssemblyGraph::createDeBruijnEdge(int node1Number, int node2Number)
{
    //Quit if any of the nodes don't exist.
    if (!m_deBruijnGraphNodes.contains(node1Number) ||
            !m_deBruijnGraphNodes.contains(node2Number) ||
            !m_deBruijnGraphNodes.contains(-node1Number) ||
            !m_deBruijnGraphNodes.contains(-node2Number))
        return;

    DeBruijnNode * node1 = m_deBruijnGraphNodes[node1Number];
    DeBruijnNode * node2 = m_deBruijnGraphNodes[node2Number];
    DeBruijnNode * negNode1 = m_deBruijnGraphNodes[-node1Number];
    DeBruijnNode * negNode2 = m_deBruijnGraphNodes[-node2Number];

    //Quit if the edge already exists
    for (size_t i = 0; i < node1->m_edges.size(); ++i)
    {
        if (node1->m_edges[i]->m_endingNode == node2)
            return;
    }

    DeBruijnEdge * forwardEdge = new DeBruijnEdge(node1, node2);
    DeBruijnEdge * backwardEdge = new DeBruijnEdge(negNode2, negNode1);

    forwardEdge->m_reverseComplement = backwardEdge;
    backwardEdge->m_reverseComplement = forwardEdge;

    m_deBruijnGraphEdges.push_back(forwardEdge);
    m_deBruijnGraphEdges.push_back(backwardEdge);

    node1->addEdge(forwardEdge);
    node2->addEdge(forwardEdge);
    negNode1->addEdge(backwardEdge);
    negNode2->addEdge(backwardEdge);
}




void AssemblyGraph::clearOgdfGraphAndResetNodes()
{
    QMapIterator<int, DeBruijnNode*> i(m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();
        i.value()->resetNode();
    }

    m_ogdfGraph->clear();
}






//http://www.code10.info/index.php?option=com_content&view=article&id=62:articledna-reverse-complement&catid=49:cat_coding_algorithms_bioinformatics&Itemid=74
QByteArray AssemblyGraph::getReverseComplement(QByteArray forwardSequence)
{
    QByteArray reverseComplement;

    for (int i = forwardSequence.length() - 1; i >= 0; --i)
    {
        char letter = forwardSequence.at(i);

        switch (letter)
        {
        case 'A': reverseComplement.append('T'); break;
        case 'T': reverseComplement.append('A'); break;
        case 'G': reverseComplement.append('C'); break;
        case 'C': reverseComplement.append('G'); break;
        case 'a': reverseComplement.append('t'); break;
        case 't': reverseComplement.append('a'); break;
        case 'g': reverseComplement.append('c'); break;
        case 'c': reverseComplement.append('g'); break;
        case 'R': reverseComplement.append('Y'); break;
        case 'Y': reverseComplement.append('R'); break;
        case 'S': reverseComplement.append('S'); break;
        case 'W': reverseComplement.append('W'); break;
        case 'K': reverseComplement.append('M'); break;
        case 'M': reverseComplement.append('K'); break;
        case 'r': reverseComplement.append('y'); break;
        case 'y': reverseComplement.append('r'); break;
        case 's': reverseComplement.append('s'); break;
        case 'w': reverseComplement.append('w'); break;
        case 'k': reverseComplement.append('m'); break;
        case 'm': reverseComplement.append('k'); break;
        case 'B': reverseComplement.append('V'); break;
        case 'D': reverseComplement.append('H'); break;
        case 'H': reverseComplement.append('D'); break;
        case 'V': reverseComplement.append('B'); break;
        case 'b': reverseComplement.append('v'); break;
        case 'd': reverseComplement.append('h'); break;
        case 'h': reverseComplement.append('d'); break;
        case 'v': reverseComplement.append('b'); break;
        case 'N': reverseComplement.append('N'); break;
        case 'n': reverseComplement.append('n'); break;
        case '.': reverseComplement.append('.'); break;
        case '-': reverseComplement.append('-'); break;
        case '?': reverseComplement.append('?'); break;
        }
    }

    return reverseComplement;
}




void AssemblyGraph::resetEdges()
{
    for (size_t i = 0; i < m_deBruijnGraphEdges.size(); ++i)
        m_deBruijnGraphEdges[i]->reset();
}


double AssemblyGraph::getMeanDeBruijnGraphCoverage(bool drawnNodesOnly)
{
    int nodeCount = 0;
    long double coverageSum = 0.0;
    long long totalLength = 0;

    QMapIterator<int, DeBruijnNode*> i(m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();
        DeBruijnNode * node = i.value();

        if (drawnNodesOnly && !node->m_drawn)
            continue;

        if (node->m_number > 0)
        {
            ++nodeCount;
            totalLength += node->m_length;
            coverageSum += node->m_length * node->m_coverage;
        }
    }

    if (totalLength == 0)
        return 0.0;
    else
        return coverageSum / totalLength;
}

double AssemblyGraph::getMaxDeBruijnGraphCoverageOfDrawnNodes()
{
    double maxCoverage = 1.0;

    QMapIterator<int, DeBruijnNode*> i(m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();

        if (i.value()->m_graphicsItemNode != 0 && i.value()->m_coverage > maxCoverage)
            maxCoverage = i.value()->m_coverage;
    }

    return maxCoverage;
}


void AssemblyGraph::resetNodeContiguityStatus()
{
    QMapIterator<int, DeBruijnNode*> i(m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();
        i.value()->m_contiguityStatus = NOT_CONTIGUOUS;
    }
}

void AssemblyGraph::resetAllNodeColours()
{
    QMapIterator<int, DeBruijnNode*> i(m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();
        if (i.value()->m_graphicsItemNode != 0)
            i.value()->m_graphicsItemNode->setNodeColour();
    }
}

void AssemblyGraph::clearAllBlastHitPointers()
{
    QMapIterator<int, DeBruijnNode*> i(m_deBruijnGraphNodes);
    while (i.hasNext())
    {
        i.next();
        DeBruijnNode * node = i.value();
        node->m_blastHits.clear();
    }
}