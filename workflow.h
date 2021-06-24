#ifndef WORKFLOW_H_
#define WORKFLOW_H_

#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include "argparse.h"
#include "bamsubset.h"
#include <iostream>

using namespace seqan;

//Checking parameters given by user
int parseBCSubsetParams(Parameters & params, std::unordered_set<std::string> & wlBarcodes, int argc, char const * argv[])
{
    if(parseCommandLine(params, argc, argv) != ArgumentParser::PARSE_OK)
    {
        std::cerr << "ERROR: Could not parse command line\n";
        return 1;
    }
    if (params.bcWlFileName == "")
    {
        std::cerr << "ERROR: whitelist file not specified. Please use option -f\n";
        return 1;
    }
    if (params.outBamFileName == "")
    {
        std::cerr << "ERROR: Output file not specified. Please use option -o\n";
        return 1;
    }
    if(!readWhitelist(wlBarcodes, params.bcWlFileName))
    {
        std::cerr << "ERROR: Could not read " << params.bcWlFileName << "\n";
        return 1;
    }
    return 0;
}

// Generate BAM subset of reads with whitelisted barcodes
int bamSubset(int argc, char const * argv[])
{
    Stats stats;
    Parameters params;
    int res = checkParser(parseCommandLine(params, argc, argv));
    if (res >= 0)
        return res;

    std::unordered_set<std::string> wlBarcodes;

    readWhitelist(wlBarcodes, params.bcWlFileName);

    // Open BamFileIn for reading
    BamFileIn inFile;
    if (!open(inFile, toCString(params.bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << params.bamFileName << " for reading.\n";
        return 1;
    }

    // Open output file BamFileOut
    std::fstream outStream;
    outStream.open(toCString(params.outBamFileName), std::ios::out);
    if (!outStream.good())
        SEQAN_THROW(FileOpenError(toCString(params.outBamFileName)));
    
    BamFileOut bamFileOut(context(inFile), outStream, Bam());

    // Access header
    BamHeader header;
    readHeader(header, inFile);

    // Write header
    processHeader(header, bamFileOut, argv);

    processBam(inFile, bamFileOut, wlBarcodes, params.bctag, params.trimming, stats);

    std::cout << "[bcsubset] Output file has been written to \'" << params.outBamFileName << "\'." << std::endl; 

    stats.report();

    return 0;
}

#endif /* WORKFLOW_H_ */
