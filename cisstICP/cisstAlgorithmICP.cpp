// ****************************************************************************
//
//    Copyright (c) 2014, Seth Billings, Russell Taylor, Johns Hopkins University
//    All rights reserved.
//
//    Redistribution and use in source and binary forms, with or without
//    modification, are permitted provided that the following conditions are
//    met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//    3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ****************************************************************************

#include "cisstAlgorithmICP.h"


#ifdef ValidateCovTreeSearch
std::ofstream validFS("../ICP_TestData/LastRun/debugCovTreeSearchValidation.txt");
vct3 validPoint;
vct3 validNorm;
int  validDatum;
double validDist;
double validAng;
double validError;
double searchError;
unsigned int numValidDatums;
unsigned int numInvalidDatums;
double validPercent;
double doubleEps = 1e-16;
int validIter = 0;
#endif

#ifdef SaveMatchesToFile
std::string saveMatchesDir("../ICP_TestData/LastRun");
int saveMatchesIter = 0;
#endif

namespace {

    // declerations
    void ICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData);

    // default callback
    void ICPCallback_PrintIteration(cisstICP::CallbackArg &arg, void *userData)
    {
        vctRodRot3 dR(arg.dF.Rotation());
        std::stringstream ss;
        cisstAlgorithmICP *pThis = (cisstAlgorithmICP*)userData;

        std::string fstring("i=%u E=%.1f tolE=%.3f  t=%.3f NNodes=%u/%u/%u NOut=%u");
        ss << cmnPrintf(fstring.c_str())
            // \tt=%.3f\tNNodes=%u/%u/%u
            // (RMS/Res)=%.4f/%.4f
            // (dAng/dPos)= %.2f/%.2f )
           << arg.iter
           << arg.E
           << arg.tolE
           << arg.time
           << pThis->maxNodesSearched << pThis->avgNodesSearched << pThis->minNodesSearched
           << arg.nOutliers
            //<< dR.Norm()*180/cmnPI << arg.dF.Translation().Norm()
            ;
#ifdef ValidateCovTreeSearch
        fstring = "  vld=%.2f";
        ss << cmnPrintf(fstring.c_str())
           << validPercent;
#endif

        std::cout << ss.str() << std::endl;
    }

} // namespac anonymous

std::vector<cisstICP::Callback> cisstAlgorithmICP::ICP_GetIterationCallbacks()
{
    std::vector<cisstICP::Callback> callbacks;
    cisstICP::Callback defaultICPCallback(ICPCallback_PrintIteration, this);
    callbacks.push_back(defaultICPCallback);
    return callbacks;
}

// constructor
cisstAlgorithmICP::cisstAlgorithmICP(cisstCovTreeBase *pTree, vctDynamicVector<vct3> &samplePts)
    : pTree(pTree)
{
    avgNodesSearched = 0;
    SetSamples(samplePts);
}

void cisstAlgorithmICP::SetSamples(vctDynamicVector<vct3> &argSamplePts)
{
    // copy sample points
    samplePts = argSamplePts;

    // initialize variables dependent on sample size
    nSamples = samplePts.size();

    samplePtsXfmd.SetSize(nSamples);
    matchPts.SetSize(nSamples);
    matchDatums.SetSize(nSamples);
    matchErrors.SetSize(nSamples);

    residuals_PostMatch.SetSize(nSamples);
    sqrDist_PostMatch.SetSize(nSamples);
    dist_PostMatch.SetSize(nSamples);

    residuals_PostRegister.SetSize(nSamples);
    sqrDist_PostRegister.SetSize(nSamples);
    dist_PostRegister.SetSize(nSamples);
}

void cisstAlgorithmICP::ICP_InitializeParameters(vctFrm3 &FGuess)
{
    // set starting sample positions
    UpdateSamplePositions(FGuess);

    // initialize matches with accelerated approximate search
    for (unsigned int i = 0; i < nSamples; i++)
    {
        matchDatums[i] = pTree->FastInitializeProximalDatum(
            samplePtsXfmd[i], matchPts[i]);
    }

    //// initialize matches to any model point
    ////   i.e. we don't know the closest match => set it to anything valid
    //for (unsigned int s = 0; s < nSamples; s++)
    //{
    //  matchDatums.Element(s) = 0;
    //  matchPts.Element(s) = pTree->DatumSortPoint(0);
    //}

    nOutliers = 0;
    Freg = FGuess;
}

void cisstAlgorithmICP::ICP_UpdateParameters_PostMatch(int index)
{
    //ComputeErrors_PostMatch();
}

void cisstAlgorithmICP::ICP_UpdateParameters_PostRegister(vctFrm3 &Freg,
                                                          int index)
{
    UpdateSamplePositions(Freg, index);

    //ComputeErrors_PostRegister();
}


/*
  ParameterizedTest_PointCloud_SurfaceNoise takes in params from main.cpp
  It then calls GenerateSamplePointSet_PointCloud_SurfaceNoies to generate the sample
  makes the target, source, and noisy target mesh

  Then it builds the cissCovTree_PointCloud ptree from noisy mesh
  should make a cisstAlgorithmICP_standalone file, maybe isolate and import
  just the files I need for it.

  need a fn that takes in a bunch of points (and makes a tree out of them)
  and a single point (and matches it) and rot/trans errors (and improves them
*/
void cisstAlgorithmICP::ICP_ComputeMatches(unsigned int &nodesSearched,
                                           int index)
{
    // Find the point on the model having lowest match error
    //  for each sample point

#ifdef ValidateCovTreeSearch
    numInvalidDatums = 0;
    numValidDatums = 0;
#endif

    minNodesSearched = std::numeric_limits<unsigned int>::max();
    maxNodesSearched = std::numeric_limits<unsigned int>::min();

    if (index <= 0)
        avgNodesSearched = 0;

    if (index == -1) {
        for (unsigned int s = 0; s < nSamples; s++)
        {
            ICP_MatchPoint(s, nodesSearched);
        }
    }
    else {
        ICP_MatchPoint(index, nodesSearched);
    }

    if (index == -1 || index == nSamples - 1)
        avgNodesSearched /= nSamples;

#ifdef ValidateCovTreeSearch
    validPercent = (double)numValidDatums / (double)nSamples;
    validFS << "iter " << validIter << ":  NumMatches(valid/invalid): "
            << numValidDatums << "/" << numInvalidDatums << "  valid% = "
            << validPercent << std::endl;
    validIter++;
#endif

#ifdef SaveMatchesToFile
    {
        std::stringstream ss;
        ss << saveMatchesDir + "/matches-" << saveMatchesIter << ".txt";
        std::ofstream fsCP(ss.str().c_str());

        ss.str("");
        ss << saveMatchesDir + "/samplesXfmd-" << saveMatchesIter << ".txt";
        std::ofstream fsSP(ss.str().c_str());

        unsigned int start = 0;
        unsigned int end = nSamples;
        if (index != -1) {
            start = index;
            end = start + 1;
        }
        for (unsigned int i = start; i < end; i++)
        {
            fsCP << matchPts.at(i) << std::endl;
            fsSP << samplePtsXfmd.at(i) << std::endl;
        }
        fsCP.close();
        fsSP.close();
        saveMatchesIter++;
    }
#endif

}

void cisstAlgorithmICP::ICP_MatchPoint(unsigned int s, unsigned int &nodesSearched)
{

    // inform algorithm beginning new match
    SamplePreMatch(s);

    // Find best match for this sample
    matchDatums.Element(s) = pTree->FindClosestDatum(
        samplePtsXfmd.Element(s), // vct3 vector
        matchPts.Element(s), // best closest point estimate so far, vct3 vector
        matchDatums.Element(s), // int prevDatum
        matchErrors.Element(s), // double matchError
        nodesSearched); //unsigned int &numNodesSearched
    avgNodesSearched += nodesSearched;
    minNodesSearched = (nodesSearched < minNodesSearched) ? nodesSearched : minNodesSearched;
    maxNodesSearched = (nodesSearched > maxNodesSearched) ? nodesSearched : maxNodesSearched;

    // these just check every datum in the tree to make sure the tree isn't broken
#ifdef ValidateCovTreeSearch
#ifdef ValidateByEuclideanDist
    validDatum = pTree->ValidateClosestDatum_ByEuclideanDist(samplePtsXfmd.Element(s), validPoint);
#else
    validDatum = pTree->ValidateClosestDatum(samplePtsXfmd.Element(s), validPoint);
#endif
    if (validDatum != matchDatums.Element(s))
    {
        // It is possible to have different datums for same point if the match
        //  lies on a datum edge; if this is the case, then the search did not
        //  actually fail
        searchError = (validPoint - matchPts.Element(s)).NormSquare();
        // Note: cannot compare two double values for exact equality due to
        //       inexact representation of decimal values in binary arithmetic
        //       => use an epsilon value for comparison

        // if the search match error was worse than the match error found by
        // checking every datum in the tree
        if (searchError > doubleEps)
        {
            double ResidualDistance = (matchPts.Element(s) - samplePtsXfmd.Element(s)).Norm();
            validDist = (validPoint - samplePtsXfmd.Element(s)).Norm();
            numInvalidDatums++;
            vct3 tmp1;

            // find point on matchDatums.Element(s) with lowest match error
            // returns match error and tmp1 set to lowest error point

            searchError = pTree->pAlgorithm->FindClosestPointOnDatum(
                samplePtsXfmd.Element(s), tmp1, matchDatums.Element(s));
            // this is the error with the point found by traversing the entire tree (valid datum)
            validError = pTree->pAlgorithm->FindClosestPointOnDatum(
                samplePtsXfmd.Element(s), tmp1, validDatum);

            validFS << "Match Errors = " << searchError << "/" << validError
                    << "\t\tdPos = " << ResidualDistance << "/" << validDist << std::endl;

            validFS << " XfmSamplePoint = [" << samplePtsXfmd.Element(s) << "]" << std::endl;
            validFS << " SearchPoint = [" << matchPts.Element(s) << "]" << std::endl;
            validFS << " SearchDatum = " << matchDatums.Element(s) << std::endl;
            validFS << " ValidPoint =  [" << validPoint << "]" << std::endl;
            validFS << " ValidDatum = " << validDatum << std::endl;
            validFS << " SampleIndex = " << s << std::endl;

            cisstCovTreeNode *termNode = 0;
            pTree->FindTerminalNode(validDatum, &termNode);
            if (!termNode)
            {
                std::cout << "ERROR: did not find terminal node for datum: " << validDatum << std::endl;
                assert(0);
            }
            validFS << " Valid Terminal Node:" << std::endl;
            validFS << "   MinCorner: " << termNode->Bounds.MinCorner << std::endl;
            validFS << "   MaxCorner: " << termNode->Bounds.MaxCorner << std::endl;
            validFS << "   NData: " << termNode->NData << std::endl;
            validFS << "Fnode = [" << std::endl << termNode->F << "]" << std::endl;

            // crash program
            std::cout << "Invalid Node Search; Terminating Program" << std::endl;
            termNode = NULL;
            termNode->NData = 0;
        }
        else
        {
            numValidDatums++;
        }
    }
    else
    {
        numValidDatums++;
    }
#endif

    // inform algorithm that match completed
    // Does this function even do anything?
    SamplePostMatch(s);
}


unsigned int cisstAlgorithmICP::ICP_FilterMatches(int index)
{
    return 0;
}

void cisstAlgorithmICP::UpdateSamplePositions(const vctFrm3 &F,
                                              int index)
{
    unsigned int start = 0;
    unsigned int end = nSamples;
    if (index != -1) {
        start = index;
        end = start + 1;
    }
    for (unsigned int s = start; s < end; s++)
    {
        samplePtsXfmd.Element(s) = F * samplePts.Element(s);
    }
}

void cisstAlgorithmICP::ComputeErrors_PostMatch(int index)
{
    //matchErrorAvg_PostMatch = 0.0;
    if (index <= 0) {
        matchDistAvg_PostMatch = 0.0;
        sumSqrDist_PostMatch = 0.0;
    }
    unsigned int start = 0;
    unsigned int end = nSamples;
    if (index != -1) {
        start = index;
        end = index + 1;
    }

    for (unsigned int s = start; s < end; s++)
    {
        residuals_PostMatch.Element(s) = samplePtsXfmd.Element(s) - matchPts.Element(s);
        sqrDist_PostMatch.Element(s) = residuals_PostMatch.Element(s).NormSquare();
        dist_PostMatch.Element(s) = sqrt(sqrDist_PostMatch.Element(s));

        sumSqrDist_PostMatch += sqrDist_PostMatch.Element(s);
        matchDistAvg_PostMatch += dist_PostMatch.Element(s);
        //matchErrorAvg_PostMatch += matchErrors.Element(s);
    }

    //matchErrorAvg_PostMatch /= nSamples;
    if (index == -1 || index == nSamples - 1)
        matchDistAvg_PostMatch /= nSamples;
}

void cisstAlgorithmICP::ComputeErrors_PostRegister(int index)
{
    //matchErrorAvg_PostRegister = 0.0;
    if (index <= 0) {
        matchDistAvg_PostRegister = 0.0;
        sumSqrDist_PostRegister = 0.0;
    }
    unsigned int start = 0;
    unsigned int end = nSamples;
    if (index != -1) {
        start = index;
        end = start + 1;
    }

    for (unsigned int s = start; s < end; s++)
    {
        residuals_PostRegister.Element(s) = samplePtsXfmd.Element(s) - matchPts.Element(s);
        sqrDist_PostRegister.Element(s) = residuals_PostRegister.Element(s).NormSquare();
        dist_PostRegister.Element(s) = sqrt(sqrDist_PostRegister.Element(s));

        sumSqrDist_PostRegister += sqrDist_PostRegister.Element(s);
        matchDistAvg_PostRegister += dist_PostRegister.Element(s);
        //matchErrorAvg_PostRegister += matchErrors.Element(s);
    }

    //matchErrorAvg_PostRegister /= nSamples;
    if (index == -1 || index == nSamples - 1)
        matchDistAvg_PostRegister /= nSamples;
}
