#ifndef WORKFLOW_H_
#define WORKFLOW_H_

#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include "argparse.h"
#include "bamsubset.h"
#include <iostream>
#include <iterator>
#include <map>

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

// Read whitelist file in

// Generate bam subset of reads with whitelisted barcodes

int bamSubset(int argc, char const * argv[])
{
    Stats stats;
    Parameters params;
    parseCommandLine(params, argc, argv);
    std::unordered_set<std::string> wlBarcodes;

    readWhitelist(wlBarcodes, params.bcWlFileName);

    // Open BamFileIn for reading.
    BamFileIn inFile;
    if (!open(inFile, toCString(params.bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << params.bamFileName << " for reading.\n";
        return 1;
    }

    // Open output file, BamFileOut accepts also an ostream and a format tag.
    
    std::fstream outStream;
    outStream.open(toCString(params.outBamFileName), std::ios::out);
    if (!outStream.good())
        SEQAN_THROW(FileOpenError(toCString(params.outBamFileName)));

    //std::cout << "Opened Stream" <<std::endl;
    
    BamFileOut bamFileOut(context(inFile), outStream, Bam()); //output file for 

    // Access header
    BamHeader header;
    readHeader(header, inFile);

    // Write header
    processHeader(header, bamFileOut, argv);


    //std::cout << "Processed Header" <<std::endl;

    processBam(inFile, bamFileOut, wlBarcodes, stats);

    stats.report();

    std::cout << "\nOutput file has been written to \'" << params.outBamFileName << "\'." << std::endl; 

    return 0;
}

#endif /* WORKFLOW_H_ */
