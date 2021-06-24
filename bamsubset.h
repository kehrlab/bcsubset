#ifndef BAMSUBSET_H_
#define BAMSUBSET_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <iostream>
#include <unordered_set>

using namespace seqan;

struct Stats
{
    unsigned filteredReads;
    unsigned passedReads;

    Stats(): filteredReads(0), passedReads(0){}

    inline void report()
    {
        std::cout << "\nSUMMARY" << std::endl;
        std::cout << "Total records:\t\t" << (filteredReads + passedReads) << std::endl;
        std::cout << "Filtered records:\t" << filteredReads << "\t(" << static_cast<double>(filteredReads)/(filteredReads + passedReads)*100 << "%)" 
                    << "\nPassed records:\t\t" << passedReads << "\t(" << static_cast<double>(passedReads)/(filteredReads + passedReads)*100 << "%)" << std::endl;
    }
};  

//Read a text file containing whitelisted barcodes, put into set of strings
//Return false if text file can not be opened
//Return true if success
bool readWhitelist(std::unordered_set<std::string> & wlBarcodes, const CharString & bcWlFileName)
{
    // Read whitelisted barcodes file
    std::ifstream wlIn(toCString(bcWlFileName));

    if (!wlIn.is_open())
    {
        std::cerr << "ERROR: Could not open " << bcWlFileName << " for reading.\n";
        return false;
    }

    std::string barcode;
    
    while (wlIn >> barcode)
    {
        // Read the file line by line, save each barcode
        wlBarcodes.emplace(barcode);
    }

    std::cout << "\n[bcsubset] Loaded " << wlBarcodes.size() << " barcodes from \'" << bcWlFileName << "\'." << std::endl;

    return !empty(wlBarcodes);
}

// Process BAM header, add @PG line
inline void processHeader(BamHeader & header, BamFileOut & bamFileOut, char const ** argv)

{
    // Modify @PG
    BamHeaderRecord hrecord;
    hrecord.type = BamHeaderRecordType::BAM_HEADER_PROGRAM;

    // Set tags
    appendValue(hrecord.tags, Pair<CharString>("ID", "bcsubset"));
    appendValue(hrecord.tags, Pair<CharString>("VN", VERSION));

    CharString clstring;
    append(clstring, "bcsubset");

    for (unsigned j=1; j<length(argv);++j)
    {
        appendValue(clstring, ' ');
        append(clstring, argv[j]);
    }

    appendValue(hrecord.tags, Pair<CharString>("CL", clstring));

    appendValue(header, hrecord);

    writeHeader(bamFileOut, header);
}

//Trim the last n characters of barcode
inline void trimBarcode(CharString & barcode, const unsigned toTrim)
{
    resize(barcode, length(barcode) - toTrim);
}

// Get barcode from tags in BAM records
inline bool getBarcodeFromTags(std::string & barcode, const BamAlignmentRecord & record, const CharString & bctag, const unsigned toTrim)
{
    unsigned idx = 0;

    BamTagsDict tagsDict(record.tags);
    bool keyFound = findTagKey(idx, tagsDict, bctag);
    if (keyFound)
    {
        CharString tag;
        if (!extractTagValue(tag, tagsDict, idx))
        {
        std::cerr << "WARNING: There was an error extracting barcode from tag " << bctag << " of record: " << record.qName << "\n";
        return false;
        }
        else
        {
        trimBarcode(tag, toTrim);
        move(barcode, tag);
        return true;
        }
    }
    else
    {
        return false;
  }
}

// Check if a BAM record contains a whitelisted barcode
inline bool isGoodRecord(const BamAlignmentRecord & record, const std::unordered_set<std::string> & wlBarcodes, const CharString & bctag, const unsigned toTrim)
{
    std::string readBC;
    if(!getBarcodeFromTags(readBC, record, bctag, toTrim))
        return false;

    if (wlBarcodes.find(readBC) == wlBarcodes.end())
        return false;
    else
        return true;
}

// Process input BAM file to find records matching the whitelisted barcodes and write them to output BAM file
inline void processBam(BamFileIn & inFile, BamFileOut & bamFileOut, const std::unordered_set<std::string> & wlBarcodes, const CharString & bctag, const unsigned toTrim, Stats & stats)
{
    while (!atEnd(inFile))
    {
        BamAlignmentRecord record;
        readRecord(record, inFile);

        if(isGoodRecord(record, wlBarcodes, bctag, toTrim))
        {
            writeRecord(bamFileOut, record);
            ++stats.passedReads;
        }
        else
        {
            ++stats.filteredReads;
        }
    }
}

#endif /* BAMSUBSET_H_ */