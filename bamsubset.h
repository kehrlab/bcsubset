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
       std::cout << "Total records:\t" << (filteredReads + passedReads) << std::endl;
       std::cout << "Filtered records:\t" << filteredReads << "\nPassed records:\t" << passedReads << std::endl;
   }

};

//Read a text file containing whitelisted barcodes, put into set of strings
//Return false if text file can not be opened
//Return true if success

bool readWhitelist(std::unordered_set<std::string> & wlBarcodes, const CharString & bcWlFileName)
{
    // Read whitelist File.
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

    std::cout << "Loaded " << wlBarcodes.size() << " barcodes from \'" << bcWlFileName << "\'." << std::endl;
    std::cout << std::endl;

    return !empty(wlBarcodes);
}


// Process bam header, add @PG line
inline void processHeader(BamHeader & header, BamFileOut & bamFileOut, char const ** argv)

{
  // Modify @PG
  BamHeaderRecord hrecord;
  hrecord.type = BamHeaderRecordType::BAM_HEADER_PROGRAM;
  
  // Set tags
  appendValue(hrecord.tags, Pair<CharString>("ID", "bcsubset"));
  appendValue(hrecord.tags, Pair<CharString>("VN", "0.0.1"));

  CharString clstring;
  append(clstring, "bcsubset ");

  for (unsigned j=1; j<length(argv);++j)
  {
    appendValue(clstring, ' ');
    append(clstring, argv[j]);
  }

  //appendValue(hrecord.tags, Pair<CharString>("CL", "bcsubset test"));
  appendValue(hrecord.tags, Pair<CharString>("CL", clstring));
  
  appendValue(header, hrecord);

  writeHeader(bamFileOut, header);
}

// Get barcode from tags in bam records, CR tag
inline void getBarcodeFromTags(std::string & barcode, const BamAlignmentRecord & record)
{
  unsigned idx = 0;

  BamTagsDict tagsDict(record.tags);
  bool keyFound = findTagKey(idx, tagsDict, "CR");
  if (keyFound)
  {
    CharString tag;
    if (!extractTagValue(tag, tagsDict, idx))
    {
      std::cerr << "ERROR: There was an error extracting CR from tags of record: " << record.qName << "\n";
    }
    move(barcode, tag);
  }
  else
  {
    std::cerr << "ERROR: CR not found in tags of record: " << record.qName << "\n";
  }
} 

inline bool isGoodRecord(const BamAlignmentRecord & record, const std::unordered_set<std::string> & wlBarcodes)
{
    std::string readBC;
    getBarcodeFromTags(readBC, record);

    if (wlBarcodes.find(readBC) == wlBarcodes.end())
        return false;
    else
        return true;
}

// Find reads that come from a whitelisted barcode, go over every record to check barcode, if it matches the whitelisted barcode, write to subset output file
//Return false if no read has the barcode
//Return true on success

inline void processBam(BamFileIn & inFile, BamFileOut & bamFileOut, const std::unordered_set<std::string> & wlBarcodes, Stats & stats)
{
    while (!atEnd(inFile))
    {
        BamAlignmentRecord record;
        readRecord(record, inFile);

        if(isGoodRecord(record, wlBarcodes))
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