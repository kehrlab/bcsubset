#ifndef BAMSUBSET_H_
#define BAMSUBSET_H_

#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include "argparse.h"
#include <iostream>
#include <iterator>
#include <map>

using namespace seqan;
using namespace std;

struct ConsensusRecord
{
    std::string fwdSeq;
    std::string revSeq;
    unsigned fwdBeginPos;
    unsigned revBeginPos;
    CharString barcode;
    unsigned seqNum;
    String<CharString> idList;
    int rID;
    int fwdAQual; //forward consensus alignment quality, MAPQ
    int revAQual; //reverse consensus alignment quality
    CharString fwdQString;  // quality string
    CharString revQString;  // quality string

    ConsensusRecord():fwdSeq(""), revSeq(""), fwdBeginPos(0), revBeginPos(0), barcode(""), seqNum(0), rID(0), fwdAQual(0), revAQual(0), fwdQString(""), revQString(""){}
};

// Jumping to a region between a start and end position
// Return -1 on error
// Return 1 on success
// Return 0 if no alignments on region 

int jumping (BamFileIn & inFile, const BamIndex<Bai> & baiIndex, int rID, int beginPos, int endPos)
{
    // 1-based to 0-based.
    beginPos -= 1;
    endPos -= 1;

    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inFile, hasAlignments, rID, beginPos, endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << beginPos << ":" << endPos << "\n";
        return -1;
    }
    if (!hasAlignments)
        return 0;  // No alignments here.

    return 1;
}

//Overloading function to use with bedRecord

int jumping(BamFileIn & inFile, const BamIndex<Bai> & baiIndex, const BedRecord<Bed3> & bedRecord)
{
    return jumping(inFile, baiIndex, bedRecord.rID, bedRecord.beginPos, bedRecord.endPos);
}

//Sequentially advance to a position in the bam file
//Record will be the first record within the region
//Return false if there is nothing there
//Return true on success

bool gotoPos(BamAlignmentRecord & record, BamFileIn & inFile, const int & rID, int beginPos, int endPos)
{
    --beginPos;
    --endPos;

    while (!atEnd(inFile))
    {
        readRecord(record, inFile);

        //std::cout<< "Record rID is "<< record.rID <<"\trecord.beginPos "<< record.beginPos <<std::endl;

        // If we are on the next reference or at the end already then we stop.
        if (record.rID == -1 || record.rID > rID || record.beginPos >= endPos)
            return false;
        // If we are left of the selected position then we skip this record.
        if (record.beginPos < beginPos)
            continue;
        return true;
    }
    return false;
}

//Overloading function to use with bedRecord

bool gotoPos(BamAlignmentRecord & record, BamFileIn & inFile, const BedRecord<Bed3> & bedRecord)
{
    return gotoPos(record, inFile, bedRecord.rID, bedRecord.beginPos, bedRecord.endPos);
}


// Check if a bam record is a good forward read based on sam flags

bool isGoodForward(const BamAlignmentRecord & record)
{
    if ((record.flag & 33) != 33)              // If read is paired and mate is in reverse strand, first in pair
        return false;
    //if ((record.flag & 3868) != 0)             // 3868: read unmapped, mate unmapped, not primary alignment, read fails quality checks, PCR duplicate, supplementary alignment or reverse strand 
    if ((record.flag & 2844) != 0)      // 2844: read unmapped, mate unmapped, read reverse strand, not primary alignment, read fails quality checks, supplementary alignment
        return false;
    if (record.pNext > record.beginPos + 10000) // Mate not to far off
        return false;
    return true;
}

// Check if a bam record is a good reverse read based on sam flags

bool isGoodReverse(const BamAlignmentRecord & record)
{
    if((record.flag & 17) != 17)              // If read is paired, in reverse strand and second in pair
        return false;
    //if ((record.flag & 3852) != 0)             // 3852: read unmapped, mate unmapped, not primary alignment, read fails quality checks, PCR duplicate, supplementary alignment
    if ((record.flag & 2828) != 0)          // 2828: read unmapped, mate unmapped, not primary aligmnent, read fails quality checks, supplementary alignment
        return false;
    return true;
}

// Process good reverse: find its forward on the map, if found, delete it and add the pair to goodPairs string
// Return true if forward is found in map, removed from it and added to goodPairs along with the reverse read
// Return false is forward read is not found in map, means good reverse read has no good forward read

bool processGoodReverse(String<Pair<BamAlignmentRecord> > & goodPairs, map<CharString, BamAlignmentRecord> & recordMap, const BamAlignmentRecord & record)
{
    map<CharString, BamAlignmentRecord>::iterator fwdIt;
    fwdIt = recordMap.find(record.qName);
    if (fwdIt == recordMap.end())
    {
        return false;
    }
    //std::cout<< "Found mate of qName:" << record.qName <<std::endl;
    seqan::appendValue(goodPairs, Pair<BamAlignmentRecord>(fwdIt->second, record)); //add good pair to goodPairs string of pairs
    recordMap.erase(fwdIt);

    return true;
}

// Convert sequence string from BamAlignmentRecord (Iupac string) into a std::string to be used by SPOA

template<typename TSource>
void convertString(std::string & target, const TSource & source) //TODO: find more efficient conversion from Iupac to std string
{
    std::stringstream stream;
    stream << source;
    target = stream.str();
}

// Create a subset of pairs of reads with the same start and end position
// Return false if we have finished processing all the pairs
// Return true if there are still pairs to look for subset

unsigned createSubset(std::vector<std::string> & sequencesFwd,
                      std::vector<std::string> & sequencesRev,
                      String<CharString> & qualStringsFwd,
                      String<CharString> & qualStringsRev,
                      unsigned & fwdBeginPos,
                      unsigned & revBeginPos,
                      String<CharString> & idList,
                      int & rID,
                      int & fwdAQual,
                      int & revAQual, 
                      CharString & currentBC, 
                      Iterator<const String<Pair<BamAlignmentRecord> > >::Type & start, 
                      const Iterator<const String<Pair<BamAlignmentRecord> > >::Type & ending)
{
    int currentPos = start->i1.beginPos;
    getBarcodeFromTags(currentBC, start->i1);

    fwdBeginPos = start->i1.beginPos; // get start position of forward reads, same as currentPos
    revBeginPos = start->i2.beginPos; // get start position of reverse reads
    rID = start->i1.rID;
    fwdAQual = start->i1.mapQ; // initialize mapQ value for forward reads, highest of subset will be chosen in the end
    revAQual = start->i2.mapQ; // initialize mapQ value for reverse reads, highest of subset will be chosen in the end

    std::string seqFwd; // store converted forward read seq
    std::string seqRev; // store converted reverse read seq
    
    CharString nextBC;

    seqan::appendValue(idList, start->i1.qName); // Initialize list of qNames of seqs making up the consensus
    convertString(seqFwd, start->i1.seq);
    sequencesFwd.push_back(seqFwd); //sequencesFwd, string vector of all forward reads seqs in std::string, pushing first sequence
    convertString(seqRev, start->i2.seq);
    sequencesRev.push_back(seqRev); //sequencesRev, string vector of all reverse reads seqs in std::string, pushing first sequence
    seqan::appendValue(qualStringsFwd, start->i1.qual); //Initialize string of CharStrings containing quality strings, forward
    seqan::appendValue(qualStringsRev, start->i2.qual); 
    ++start;

    for(; start != ending; ++start)
    {
        getBarcodeFromTags(nextBC, start->i1);
        if ((start->i1.beginPos == currentPos) && (nextBC == currentBC))
            {
                seqan::appendValue(idList, start->i1.qName);
                seqan::appendValue(qualStringsFwd, start->i1.qual); //Initialize string of CharStrings containing quality strings, forward
                seqan::appendValue(qualStringsRev, start->i2.qual);
                convertString(seqFwd, start->i1.seq); // converting and pushing from second forward read sequence on until beginPos of forward read in start is different
                sequencesFwd.push_back(seqFwd);
                convertString(seqRev, start->i2.seq); // converting and pushing from second reverse read on
                sequencesRev.push_back(seqRev);
                if (start->i1.mapQ > fwdAQual)
                    fwdAQual = start->i1.mapQ; // If higher alignment quality is found, replace fwdAQual for higher value
                if (start->i2.mapQ > revAQual)
                    revAQual = start->i2.mapQ; 
            }        
        else
            break;
    }

    return length(sequencesFwd);
}

// Custom Key Sorting function, sorting pairs of reads by starting position of the forward read and then by barcode
// Return true if the starting position and the barcode are smaller (numerically and alphabetically respectively)
// Return false otherwise

bool smallerReadPair(const Pair<BamAlignmentRecord> & left, const Pair<BamAlignmentRecord> & right)
{
    if(left.i1.beginPos < right.i1.beginPos) //sort by forwar read beginPos
        return true;
    if(left.i1.beginPos > right.i1.beginPos)
        return false;
    
    CharString leftBC; // TODO: Can be optimized to extract tags just once
    CharString rightBC;

    getBarcodeFromTags(leftBC, left.i1);
    getBarcodeFromTags(rightBC, right.i1);

    if (leftBC < rightBC)
        return true;
    else
        return false;
}

// Generate consensus sequence of good pairs with SPOA library, to be used in processGoodPairs

void makeConsensus(std::vector<std::string> & sequences, std::string & consensusSeq)
{   
    if (length(sequences) == 1)
    {
        consensusSeq = sequences[0];
        std::cout << ">Consensus LN:i:" << sequences[0].size() << std::endl
              << sequences[0] << std::endl;        // Print consensus sequence and its length
        
        return;
    }
    
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps

    spoa::Graph graph{};

    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    consensusSeq = graph.GenerateConsensus();    //Save each consensus into vector of strings consensusSeqs
    
    std::cout << ">Consensus LN:i:" << consensusSeq.size() << std::endl
              << consensusSeq << std::endl;        // Print consensus sequence and its length
    
    auto msa = graph.GenerateMultipleSequenceAlignment();   // Show the alignment, not added to bam file
    
    std::cout << "MSA:" << std::endl;

    for (const auto& it : msa) {
        std::cout << it << std::endl;
    }
    
    std::cout << "\n" << std::endl;

    return;
}

void getConsensusQuality(CharString & consQual, const std::string & consensusSeq, const std::vector<std::string> & sequences, const String<CharString> & qualSequences)
{
    SEQAN_ASSERT_EQ(length(sequences), length(qualSequences));
    
    resize(consQual, length(consensusSeq), '0', Exact());
    for (unsigned c=0; c < length(consensusSeq); ++c)
    {
        for (unsigned s=0; s < length(sequences); ++s)
        {
            SEQAN_ASSERT_EQ(length(sequences[s]), length(qualSequences[s]));
            if (c >= length(sequences[s]))
                continue;

            if (sequences[s][c] == consensusSeq[c])
            {
                if(qualSequences[s][c] > consQual[c])
                    consQual[c] = qualSequences[s][c];
            }
        }
    }
}

void convertToBamRecord(BamAlignmentRecord & fwdRecord, BamAlignmentRecord & revRecord, const ConsensusRecord & consensusRecord)
{
    // Set forward consensus bam record, sam flag not set
    fwdRecord.beginPos = consensusRecord.fwdBeginPos;
    fwdRecord.seq = consensusRecord.fwdSeq;
    fwdRecord.qName = "cons_";
    append(fwdRecord.qName, consensusRecord.idList[0]);
    fwdRecord.rID = consensusRecord.rID;
    fwdRecord.flag = 97; // paired, mate in reverse, first in pair 
    fwdRecord.pNext = consensusRecord.revBeginPos;
    fwdRecord.rNextId = consensusRecord.rID;
    fwdRecord.mapQ = consensusRecord.fwdAQual;
    fwdRecord.qual = consensusRecord.fwdQString;
    fwdRecord.tLen = consensusRecord.revBeginPos + length(consensusRecord.revSeq) - consensusRecord.fwdBeginPos + 1;
    BamTagsDict fwdTagsDict(fwdRecord.tags);
    setTagValue(fwdTagsDict, "BC", consensusRecord.barcode);
    setTagValue(fwdTagsDict, "SN", consensusRecord.seqNum);
    CharString idList;
    for (unsigned i=0; i < length(consensusRecord.idList); ++i)
    {
        append(idList, consensusRecord.idList[i]);
        if (i != length(consensusRecord.idList)-1)
            appendValue(idList, ' ');
    }
    setTagValue(fwdTagsDict, "IL", idList);

    // Set reverse consensus bam record
    revRecord.beginPos = consensusRecord.revBeginPos;
    revRecord.seq = consensusRecord.revSeq;
    revRecord.qName = fwdRecord.qName;
    revRecord.rID = consensusRecord.rID;
    revRecord.flag = 145; // paired and reverse strand, second in pair
    revRecord.pNext = consensusRecord.fwdBeginPos;
    revRecord.rNextId = consensusRecord.rID;
    revRecord.mapQ = consensusRecord.revAQual;
    revRecord.qual = consensusRecord.revQString;
    revRecord.tLen = -fwdRecord.tLen;
    BamTagsDict revTagsDict(revRecord.tags);
    setTagValue(revTagsDict, "BC", consensusRecord.barcode);
    setTagValue(revTagsDict, "SN", consensusRecord.seqNum);
    setTagValue(revTagsDict, "IL", idList);
}

// Process pairs of good reads stored in goodPairs, create a consensusRecord object
// Call sorting function smallerReadPair, call createSubset and makeConsensus
// Write consensus sequences as output to a bam file

bool processGoodPairs(String<Pair<BamAlignmentRecord>> & goodPairs, BamFileOut & bamFileOut)
{
    //std::cout<< "Processing good pairs:"<<std::endl;
    for (auto& it : goodPairs) // print good pairs before sorting
    {
        CharString barcode;
        getBarcodeFromTags(barcode, it.i1);
        //std::cout << it.i1.qName << "\t" << it.i1.beginPos << "\t" << it.i2.beginPos << "\t" << "BC: "<< barcode << "\tMAPQ fwd: " << it.i1.mapQ << "\tMAPQ rev: " << it.i2.mapQ << std::endl;
    }

    //std::cout << goodPairs[i].i1.qName << "\t" << goodPairs[i].i1.beginPos << "\t" << goodPairs[i].i2.beginPos <<std::endl;
    std::sort(begin(goodPairs), end(goodPairs), smallerReadPair);
    
    //std::cout<< "Printing sorted pairs:"<<std::endl;

    for (auto& it : goodPairs) // print good pairs after sorting
    {
        CharString barcode;
        getBarcodeFromTags(barcode, it.i1);
        //std::cout << it.i1.qName << "\t" << it.i1.beginPos << "\t" << it.i2.beginPos << "\tBC: "<< barcode << "\tMAPQ fwd: " << it.i1.mapQ << "\tMAPQ rev: " << it.i2.mapQ << std::endl;
        std::cout << it.i1.qName << "\t" << it.i1.beginPos << "\t" << it.i2.beginPos << "\tBC: "<< barcode << "\tfwdSeq: " << it.i1.seq << "\tMAPQ fwd: " << it.i1.mapQ << "\trevSeq: " << it.i2.seq << "\tMAPQ rev: " << it.i2.mapQ << std::endl;

    }

00E0<000E000A0E000000<000<060AE<000EAA0EAEEEA0EEEE0EEEA0EEEE0EEE00E0EEAEE<6E00EEAE00EEEEEEEEEEEA00E0<000E000A0E000000<000<060AE<000EAA0EAEEEA0EEEE0EEEA0EEEE0EEE00E0EEAEE<6E00EEAE00EEEEEEEEEEEA00E0<000E000A0E000000<000<060AE<000EAA0EAEEEA0EEEE0EEEA0EEEE0EEE00E0EEAEE<6E00EEAE00EEEEEEEEEEEA00E0<000E000A0E000000<000<060AE<000EAA0EAEEEA0EEEE0EEEA0EEEE0EEE00E0EEAEE<6E00EEAE00EEEEEEEEEEEA    Iterator<const String<Pair<BamAlignmentRecord> > >::Type start = begin(goodPairs);
    Iterator<const String<Pair<BamAlignmentRecord> > >::Type ending = end(goodPairs);
    std::vector<std::string> sequencesFwd;
    std::vector<std::string> sequencesRev;
    String<CharString> qualStringsFwd;
    String<CharString> qualStringsRev;

    String<ConsensusRecord> consensusRecords;

    while(start != ending)
    {
        ConsensusRecord consensusRecord;
        consensusRecord.seqNum = createSubset(sequencesFwd, sequencesRev,
                                              qualStringsFwd, qualStringsRev,
                                              consensusRecord.fwdBeginPos, 
                                              consensusRecord.revBeginPos, 
                                              consensusRecord.idList,
                                              consensusRecord.rID,
                                              consensusRecord.fwdAQual,
                                              consensusRecord.revAQual,
                                              consensusRecord.barcode, 
                                              start, ending);  
        makeConsensus(sequencesFwd, consensusRecord.fwdSeq);
        makeConsensus(sequencesRev, consensusRecord.revSeq);
        getConsensusQuality(consensusRecord.fwdQString, consensusRecord.fwdSeq, sequencesFwd, qualStringsFwd);
        getConsensusQuality(consensusRecord.revQString, consensusRecord.revSeq, sequencesRev, qualStringsRev);
        appendValue(consensusRecords, consensusRecord); 

        clear(sequencesFwd);
        clear(sequencesRev);
        clear(qualStringsFwd);
        clear(qualStringsRev);
    }
 
    for(unsigned i = 0; i < length(consensusRecords); i++)
	{
        BamAlignmentRecord fwdRecord;
        BamAlignmentRecord revRecord;
        convertToBamRecord(fwdRecord, revRecord, consensusRecords[i]);
        writeRecord(bamFileOut, fwdRecord);
        writeRecord(bamFileOut, revRecord);
    } 

    // Print all consensus sequences for a set of good pairs
    
    std::cout << "Consensus Sequences:" << std::endl;
    for(unsigned i = 0; i < length(consensusRecords); i++)
	{
		//std::cout << consensusRecords[i].fwdSeq <<  "\tBC: " << consensusRecords[i].barcode << "\tseqNum: " << consensusRecords[i].seqNum << std::endl;
        std::cout << consensusRecords[i].fwdSeq <<  "\tBC: " << consensusRecords[i].barcode << "\tseqNum: " << consensusRecords[i].seqNum << "\tfwdBeginPos:" << consensusRecords[i].fwdBeginPos << "\trevBeginPos:" << consensusRecords[i].revBeginPos << 
        "\tMAPQ fwd: " << consensusRecords[i].fwdAQual << std::endl;
        //std::cout << "idList: " << consensusRecords[i].idList <<  std::endl; // Currently prints each qName in a new line, TODO: join together with a separatorq
        std::cout << consensusRecords[i].revSeq << "\tMAPQ rev: " << consensusRecords[i].revAQual <<std::endl;
    }
    

    clear(goodPairs);
    return true;
}

// Only checking if record.beginPos > endPos
// Return true if record.beginPos > endPos
// Return false if record.beginPos outside of range or on a new chromosome

bool isStillInRegion(const BamAlignmentRecord & record, const BedRecord<Bed3> & bedRecord)
{
    if (record.rID != bedRecord.rID)         // If we are on a new chromosome (rID), rID from bedRecords is -1 for all bedrecords, 0>-1 so it breaks, needs to ne converted with getIDname too!
        return false;
    else if (record.beginPos +1 > bedRecord.endPos) //if position is outside range bam beginPos >  range endPos, works, records added to map
        return false;
    
    SEQAN_ASSERT_GEQ(record.beginPos +1, bedRecord.beginPos); //Check that record.beginPos +1 is greater than or equal to bedRecordbeginPos, if so, function true

    return true;
}

// Check reads for conditions to add and remove records to map, write out filtered out reads to file
// Find good pairs and store in a string of lists to be processed by processGoodPairs

//ProcessRegion writing out discarded reads into bam file
void processRegion(BamAlignmentRecord & record, BamFileIn & inFile, BamFileOut & bamFileOut, BamFileOut & outReadsBamFileOut, const BedRecord<Bed3> & bedRecord)
{
    // Create an empty map
    map<CharString, BamAlignmentRecord> recordMap;

    String<Pair<BamAlignmentRecord>> goodPairs;

    int currentPos = record.beginPos;

    while(!atEnd(inFile))                                                                
    {   
        if(isStillInRegion(record, bedRecord))
        {
            //If read is forward
            if(isGoodForward(record))
            {
                recordMap.insert(pair<CharString, BamAlignmentRecord>(record.qName, record));  // If record still in region, add to map
                //std::cout<< "Record added to map:" << record.qName <<std::endl;
            } 
            else if (isGoodReverse(record))
            {
                //std::cout<< "Found good reverse:" << record.qName <<std::endl;
                //processGoodReverse(goodPairs, recordMap, record); // Find its fwd read on map, if found add to goodPairs string of pairs
                if (!processGoodReverse(goodPairs, recordMap, record)) //If good fwd read not found, filter good rev out, if found, added pair to goodPairs
                {
                    //Write to filtered out reads, good reverse but no good fwd
                    //filteredOutStream << record.qName << "\t" << record.flag << "\t" << record.beginPos << "\tGood reverse read, no good forward read found.\n";
                    //std::cout<< "Good rev, no good fwd:" << record.qName <<std::endl;
                    // Access records tags, add tag OR:
                    BamTagsDict tagsDict(record.tags);
                    setTagValue(tagsDict, "OR", 0); //OR for Out read reason, 0 is Good reverse read, no good forward read found.
                    //writeRecord(bamFileOut, record);
                    writeRecord(outReadsBamFileOut, record);
                }
            }
            else
            {
                // Write to filtered out reads, in region but no good reverse and no good forward
                //filteredOutStream << record.qName << "\t" << record.flag << "\t" << record.beginPos << "\tNot a good reverse read nor a good forward.\n";
                // Access records tags, add tag OR:
                BamTagsDict tagsDict(record.tags);
                setTagValue(tagsDict, "OR", 1); //OR for Out read reason, 1 is Not a good reverse read, nor a good forward.
                //writeRecord(bamFileOut, record);
                writeRecord(outReadsBamFileOut, record);
            }
        }
        else
        {
            if (record.beginPos +1 < bedRecord.endPos + 10000) //If not in region, check next 10K
            {
                //std::cout<< "Record within 10k of end of region:" << record.qName <<std::endl;
                if (isGoodReverse(record))
                {
                    //std::cout<< "Found good reverse in extended region:" << record.qName <<std::endl;
                    if (!processGoodReverse(goodPairs, recordMap, record))
                    {
                        // write to filtered out reads, in next 10K region, good reverse but no good fwd
                        //filteredOutStream << record.qName << "\t" << record.flag << "\t" << record.beginPos << "\tGood reverse read in region+10K, no good forward read found in region.\n";
                        // Access records tags, add tag OR:
                        BamTagsDict tagsDict(record.tags);
                        setTagValue(tagsDict, "OR", 2); //OR for Out read reason, 2 is Good reverse read in region+10K, no good forward read found in region.
                        //writeRecord(bamFileOut, record);
                        writeRecord(outReadsBamFileOut, record);

                    }
                }
            }
            else
            {
                // Write read, not in region
                break;
            }
        }

        // Read record, go to the next
        readRecord(record, inFile);
        if (record.beginPos != currentPos)
        {
            //std::cout<< "Previous currentPos:" << currentPos << "\t";
            currentPos = record.beginPos;
            //std::cout<< "New currentPos:" << currentPos <<std::endl;
            processGoodPairs(goodPairs, bamFileOut);
        }
    }

    if (!empty(recordMap))
    {   
        // Write out to filtered out reads, good forward, no good reverse
        //filteredOutStream << "Map not empty, map size: " << recordMap.size() << "\n";
        for (map<CharString, BamAlignmentRecord>::iterator mapIt = recordMap.begin(); mapIt != recordMap.end(); ++mapIt) // print map content
        {
            //filteredOutStream << mapIt->first << "\t" << mapIt->second.flag << "\t" << mapIt->second.beginPos << "\tGood forward read, no good reverse read found.\n"; //recordMap contains pairs, first is qName, second is the whole record
            //writeRecord(outReadsBamFileOut, mapIt->second);
            // Access records tags, add tag OR:
            BamTagsDict tagsDict(mapIt->second.tags);
            setTagValue(tagsDict, "OR", 3); //OR for Out read reason, 2 is Good forward read, no good reverse read found.
            //writeRecord(bamFileOut, mapIt->second);
            writeRecord(outReadsBamFileOut, mapIt->second);
        }
        
        //filteredOutStream << "\n";
    }
    
    return;
}

// ProcessRegion without writing out of discarded reads
void processRegion(BamAlignmentRecord & record, BamFileIn & inFile, BamFileOut & bamFileOut, const BedRecord<Bed3> & bedRecord)
{
    // Create an empty map
    map<CharString, BamAlignmentRecord> recordMap;

    String<Pair<BamAlignmentRecord>> goodPairs;

    int currentPos = record.beginPos;

    while(!atEnd(inFile))                                                                
    {   
        if(isStillInRegion(record, bedRecord))
        {
            //If read is forward
            if(isGoodForward(record))
            {
                recordMap.insert(pair<CharString, BamAlignmentRecord>(record.qName, record));  // If record still in region, add to map
                //std::cout<< "Record added to map:" << record.qName <<std::endl;
            } 
            else if (isGoodReverse(record))
            {
                //std::cout<< "Found good reverse:" << record.qName <<std::endl;
                //processGoodReverse(goodPairs, recordMap, record); // Find its fwd read on map, if found add to goodPairs string of pairs
                processGoodReverse(goodPairs, recordMap, record); //If good fwd read not found, filter good rev out, if found, added pair to goodPairs
            }
        }
        else
        {
            if (record.beginPos +1 < bedRecord.endPos + 10000) //If not in region, check next 10K
            {
                //std::cout<< "Record within 10k of end of region:" << record.qName <<std::endl;
                if (isGoodReverse(record))
                {
                    //std::cout<< "Found good reverse in extended region:" << record.qName <<std::endl;
                    processGoodReverse(goodPairs, recordMap, record);
                }
            }
            else
            {
                // Write read, not in region
                break;
            }
        }

        // Read record, go to the next
        readRecord(record, inFile);
        if (record.beginPos != currentPos)
        {
            //std::cout<< "Previous currentPos:" << currentPos << "\t";
            currentPos = record.beginPos;
            //std::cout<< "New currentPos:" << currentPos <<std::endl;
            processGoodPairs(goodPairs, bamFileOut);
        }
    }
    
    return;
}

//Read a bed file and create a string of records
//Return false if bed file can not be opened
//Return true if success

bool readBED(String<BedRecord<Bed3> > & bedRecords, const CharString & bedFileName)
{
    // Read BED File.
    BedFileIn bedIn;
    if (!open(bedIn, toCString(bedFileName)))
    {
        std::cerr << "ERROR: Could not open " << bedFileName << " for reading.\n";
        return false;
    }

    while (!atEnd(bedIn))
    {
        // Read the file record by record.
        BedRecord<Bed3> bedRecord;

        readRecord(bedRecord, bedIn);
        appendValue(bedRecords, bedRecord);
    }
    return !empty(bedRecords);
}

int assignrIDs(String<BedRecord<Bed3> > & bedRecords, BamFileIn & inFile)
{
    for (unsigned i=0; i<length(bedRecords);++i)
    {
        // Translate from reference name to rID.
        if (!getIdByName(bedRecords[i].rID, contigNamesCache(context(inFile)), bedRecords[i].ref))
        {
            std::cerr << "ERROR: Reference sequence named " << bedRecords[i].ref << " not known.\n";
            return -1;
        }
    }
    return 1;
}

#endif /* BAMSUBSET_H_ */