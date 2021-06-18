#ifndef WORKFLOW_H_
#define WORKFLOW_H_

#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include "argparse.h"
#include "barcode.h"
#include "bamsubset.h"
#include <iostream>
#include <iterator>
#include <map>

using namespace seqan;
using namespace std;

// Processing barcodes

int bcprocess(int argc, char const * argv[])
{

    Parameters params;

    if(parseCommandLine(params, argc, argv) != ArgumentParser::PARSE_OK)
    {
        std::cerr << "ERROR: Could not parse command line\n";
        return 1;
    }

    if (params.outBamFileName == "")
    {
        std::cerr << "ERROR: Output file not specified. Please use option -o\n";
        return 1;
    }
    
    std::cout << "[bcprocess]   Processing barcodes\n" <<std::endl;

    processBarcodes(params.bamFileName, params.outBamFileName, argv, 5);
    //processBarcodes(params.bamFileName, params.outBamFileName, argv, 15); // setting length here overrides everything in barcode.h it seems 

    std::cout << "[bcprocess]   Barcodes processed\n" <<std::endl;

    return 0;
}

int parseConsensusParams(Parameters & params, String<BedRecord<Bed3> > & bedRecords, int argc, char const * argv[])
{

    if(parseCommandLine(params, argc, argv) != ArgumentParser::PARSE_OK)
    {
        std::cerr << "ERROR: Could not parse command line\n";
        return 1;
    }

    if (params.bedFileName == "")
    {
        std::cerr << "ERROR: bed file not specified. Please use option -b\n";
        return 1;
    }
    
    if (params.outBamFileName == "")
    {
        std::cerr << "ERROR: Output file not specified. Please use option -o\n";
        return 1;
    }

    if(!readBED(bedRecords, params.bedFileName))
    {
        std::cerr << "ERROR: Could not read " << params.bedFileName << "\n";
        return 1;
    }
    return 0;
}


// Make consensus sequences
// Version of function where discarded reads are written out (filteredOutStream)

int consensus(Parameters & params, String<BedRecord<Bed3> > & bedRecords, char const * argv[])
{

    // Open BamFileIn for reading.
    BamFileIn inFile;
    if (!open(inFile, toCString(params.bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << params.bamFileName << " for reading.\n";
        return 1;
    }

    // Read BAI index.
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, toCString(params.baiFileName)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << params.baiFileName << "\n";
        return 1;
    }

    // Open output file, BamFileOut accepts also an ostream and a format tag.
    
    std::fstream outStream;
    outStream.open(toCString(params.outBamFileName), std::ios::out);
    if (!outStream.good())
        SEQAN_THROW(FileOpenError(toCString(params.outBamFileName)));

    //std::cout << "Opened Stream" <<std::endl;

    // Open output file for filtered out reads, outReadsBamFileOut
    
    std::fstream filteredOutStream;
    filteredOutStream.open(toCString(params.filteredOutBamFileName), std::ios::out); //opening stream
    if (!filteredOutStream.good())
        SEQAN_THROW(FileOpenError(toCString(params.filteredOutBamFileName)));

    BamFileOut outReadsBamFileOut(context(inFile), filteredOutStream, Bam()); //output file for filtered out reads
    
    BamFileOut bamFileOut(context(inFile), outStream, Bam()); //output file for consensus

    // Access header
    BamHeader header;
    readHeader(header, inFile);

    // Write header
    processHeader(header, bamFileOut, "consensus", argv);
    processORHeader(header, outReadsBamFileOut);

    //std::cout << "Processed Header" <<std::endl;

    //std::cout<< "Jumping in file"<<std::endl;

    if(!assignrIDs(bedRecords, inFile))
    {
        std::cerr << "ERROR: Could not assign rIDs.";
        return 1;
    }

    //Moving along string of bed records
    //Applying jumping and gotoPos to each record in the string bedRecords
    for (unsigned i=0; i<length(bedRecords);++i)
    {
        if (jumping(inFile, baiIndex, bedRecords[i]) == -1)
	        return 1; 
    
        //std::cout<< "Jumping succesful"<<std::endl;

        // Seek linearly to the selected position.
        BamAlignmentRecord record;

        // Jump to first record in range
        if (!gotoPos(record, inFile, bedRecords[i]))
            return 1;
        /*
        std::cout<< "i value: "<< i <<std::endl;
        std::cout<< "bedRecord[i].rID: " << bedRecords[i].rID << std::endl;
        std::cout<< "record.qName: " << record.qName << std::endl;
        std::cout<< "record.rID: " << record.rID << std::endl;
        std::cout<< "bedRecord[i].beginPos: " << bedRecords[i].beginPos << std::endl;
        std::cout<< "bedRecord[i].endPos: " << bedRecords[i].endPos << std::endl;
        std::cout<< "Moved to bam record with begin pos: " << record.beginPos << std::endl;
        */

        //Logging for each region
        std::cout<< "[consensus]    Processing bed region: " << bedRecords[i].rID << "\tbeginPos: " <<bedRecords[i].beginPos << "\tendPos: " << bedRecords[i].endPos << "\n";

        //filteredOutStream << "bedRecord[i].rID: " << bedRecords[i].rID << "\tbeginPos: " << bedRecords[i].beginPos << "\tendPos: " << bedRecords[i].endPos << "\n\n";
        
        // Check if reads are a good forward then add to map, otherwise check if good reverse, find its forward on map, delete it and add good pair to list of strings
        //processRegion(record, inFile, bamFileOut,bedRecords[i]);
        //processRegion(record, inFile, bamFileOut, filteredOutStream, bedRecords[i]);
        processRegion(record, inFile, bamFileOut, outReadsBamFileOut, bedRecords[i]);
    }

    return 0;
}

struct Logging
{};

//Overload of consensus function to not write discarded reads, (logging added as arg)
int consensus(Parameters & params, String<BedRecord<Bed3> > & bedRecords, char const * argv[], Logging logging)
{

    (void) logging;
    
    // Open BamFileIn for reading.
    BamFileIn inFile;
    if (!open(inFile, toCString(params.bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << params.bamFileName << " for reading.\n";
        return 1;
    }

    // Read BAI index.
    BamIndex<Bai> baiIndex;
    if (!open(baiIndex, toCString(params.baiFileName)))
    {
        std::cerr << "ERROR: Could not read BAI index file " << params.baiFileName << "\n";
        return 1;
    }

    // Open output file, BamFileOut accepts also an ostream and a format tag.
    
    std::fstream outStream;
    outStream.open(toCString(params.outBamFileName), std::ios::out);
    if (!outStream.good())
        SEQAN_THROW(FileOpenError(toCString(params.outBamFileName)));

    //std::cout << "Opened Stream" <<std::endl;
    
    BamFileOut bamFileOut(context(inFile), outStream, Bam()); //output file for consensus

    // Access header
    BamHeader header;
    readHeader(header, inFile);

    // Write header
    processHeader(header, bamFileOut, "consensus", argv);

    //std::cout << "Processed Header" <<std::endl;

    //std::cout<< "Jumping in file"<<std::endl;

    if(!assignrIDs(bedRecords, inFile))
    {
        std::cerr << "ERROR: Could not assign rIDs.";
        return 1;
    }

    //Moving along string of bed records
    //Applying jumping and gotoPos to each record in the string bedRecords
    for (unsigned i=0; i<length(bedRecords);++i)
    {
        if (jumping(inFile, baiIndex, bedRecords[i]) == -1)
	        return 1; 
    
        //std::cout<< "Jumping succesful"<<std::endl;

        // Seek linearly to the selected position.
        BamAlignmentRecord record;

        // Jump to first record in range
        if (!gotoPos(record, inFile, bedRecords[i]))
            return 1;
        /*
        std::cout<< "i value: "<< i <<std::endl;
        std::cout<< "bedRecord[i].rID: " << bedRecords[i].rID << std::endl;
        std::cout<< "record.qName: " << record.qName << std::endl;
        std::cout<< "record.rID: " << record.rID << std::endl;
        std::cout<< "bedRecord[i].beginPos: " << bedRecords[i].beginPos << std::endl;
        std::cout<< "bedRecord[i].endPos: " << bedRecords[i].endPos << std::endl;
        std::cout<< "Moved to bam record with begin pos: " << record.beginPos << std::endl;
        */

        //Logging for each region
        std::cout<< "[consensus]    Processing bed region: " << bedRecords[i].rID << "\tbeginPos: " <<bedRecords[i].beginPos << "\tendPos: " << bedRecords[i].endPos << "\n";

        //filteredOutStream << "bedRecord[i].rID: " << bedRecords[i].rID << "\tbeginPos: " << bedRecords[i].beginPos << "\tendPos: " << bedRecords[i].endPos << "\n\n";
        
        // Check if reads are a good forward then add to map, otherwise check if good reverse, find its forward on map, delete it and add good pair to list of strings
        //processRegion(record, inFile, bamFileOut,bedRecords[i]);
        //processRegion(record, inFile, bamFileOut, filteredOutStream, bedRecords[i]);
        processRegion(record, inFile, bamFileOut, bedRecords[i]);
    }

    return 0;
}

#endif /* WORKFLOW_H_ */