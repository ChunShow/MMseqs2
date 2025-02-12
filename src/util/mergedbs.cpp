// #include "IndexReader.h"
// #include "DBReader.h"
// #include "DBWriter.h"
// #include "Debug.h"
// #include "Parameters.h"
// #include "Util.h"

// int mergedbs(int argc, const char **argv, const Command& command) {
//     Parameters& par = Parameters::getInstance();
//     par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

//     if (par.filenames.size() <= 2) {
//         Debug(Debug::ERROR) << "Need at least two databases for merging\n";
//         EXIT(EXIT_FAILURE);
//     }

//     const std::vector<std::string> prefices = Util::split(par.mergePrefixes, ",");

//     const int preloadMode = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) ? IndexReader::PRELOAD_INDEX : 0;
//     IndexReader qDbr(par.db1, 1, IndexReader::SEQUENCES, preloadMode, DBReader<unsigned int>::USE_INDEX);

//     // skip par.db{1,2}
//     const size_t fileCount = par.filenames.size() - 2;
//     DBReader<unsigned int> **filesToMerge = new DBReader<unsigned int>*[fileCount];
//     for (size_t i = 0; i < fileCount; i++) {
//         std::string indexName = par.filenames[i + 2] + ".index";
//         filesToMerge[i] = new DBReader<unsigned int>(par.filenames[i + 2].c_str(), indexName.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
//         filesToMerge[i]->open(DBReader<unsigned int>::NOSORT);
//     }

//     DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed, filesToMerge[0]->getDbtype());
//     writer.open();

//     Debug(Debug::INFO) << "Merging the results to " << par.db2.c_str() << "\n";
//     Debug::Progress progress(qDbr.sequenceReader->getSize());
//     for (size_t id = 0; id < qDbr.sequenceReader->getSize(); id++) {
//         progress.updateProgress();
//         unsigned int key = qDbr.sequenceReader->getDbKey(id);
//         // get all data for the id from all files
//         writer.writeStart(0);
//         for (size_t i = 0; i < fileCount; i++) {
//             size_t entryId = filesToMerge[i]->getId(key);
//             if (entryId == UINT_MAX) {
//                 continue;
//             }
//             const char *data = filesToMerge[i]->getData(entryId, 0);
//             if (data == NULL) {
//                 if (par.mergeStopEmpty == true) {
//                     break;
//                 } else {
//                     continue;
//                 }
//             }
//             if (i < prefices.size()) {
//                 writer.writeAdd(prefices[i].c_str(), prefices[i].size(), 0);
//             }
//             writer.writeAdd(data, filesToMerge[i]->getEntryLen(entryId) - 1, 0);
//         }
//         writer.writeEnd(key, 0);
//     }
//     writer.close();
//     for (size_t i = 0; i < fileCount; i++) {
//         filesToMerge[i]->close();
//         delete filesToMerge[i];
//     }
//     delete[] filesToMerge;

//     return EXIT_SUCCESS;
// }
#include "IndexReader.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Parameters.h"
#include "Util.h"
#include <set>
int mergedbs(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    if (par.filenames.size() <= 2) {
        Debug(Debug::ERROR) << "Need at least two databases for merging\n";
        EXIT(EXIT_FAILURE);
    }

    const std::vector<std::string> prefices = Util::split(par.mergePrefixes, ",");

    const int preloadMode = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) ? IndexReader::PRELOAD_INDEX : 0;
    IndexReader qDbr(par.db1, 1, IndexReader::SEQUENCES, preloadMode, DBReader<unsigned int>::USE_INDEX);

    // skip par.db{1,2}
    const size_t fileCount = par.filenames.size() - 2;
    DBReader<unsigned int> **filesToMerge = new DBReader<unsigned int>*[fileCount];
    for (size_t i = 0; i < fileCount; i++) {
        std::string indexName = par.filenames[i + 2] + ".index";
        filesToMerge[i] = new DBReader<unsigned int>(par.filenames[i + 2].c_str(), indexName.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        filesToMerge[i]->open(DBReader<unsigned int>::NOSORT);
    }

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), 1, par.compressed, filesToMerge[0]->getDbtype());
    writer.open();
    char key[255];
    char buffer[1024 + 32768 * 4];
    std::string resultOutString;
    resultOutString.reserve(1024);
    Debug(Debug::INFO) << "Merging the results to " << par.db2.c_str() << "\n";
    Debug::Progress progress(qDbr.sequenceReader->getSize());
    // writer.open();
    for (size_t id = 0; id < qDbr.sequenceReader->getSize(); id++) {
        progress.updateProgress();
        unsigned int key = qDbr.sequenceReader->getDbKey(id);
        // std::cout << "referencekey: " << key << std::endl;
        
        // get all data for the id from all files
        unsigned int repKey = key;
        resultOutString.append(SSTR(repKey));
        resultOutString.push_back('\n');
        std::set<unsigned int> memkeyset;
        // writer.writeStart(0);
        for (size_t i = 0; i < fileCount; i++) {
            size_t entryId = filesToMerge[i]->getId(key);
            // std::cout << "entryId: " << entryId << std::endl;
            if (entryId == UINT_MAX) {
                continue;
            }
            char *data = filesToMerge[i]->getData(entryId, 0);
            if (data == NULL) {
                if (par.mergeStopEmpty == true) {
                    break;
                } else {
                    continue;
                }
            }

            while (*data != '\0') {
                Util::parseKey(data, buffer);
                const unsigned int key = (unsigned int) strtoul(buffer, NULL, 10);
                // std::cout << "member key: " << key << std::endl;
                memkeyset.insert(key);
                data = Util::skipLine(data);
            }
            
            // if (i < prefices.size()) {
            //     writer.writeAdd(prefices[i].c_str(), prefices[i].size(), 0);
            // }
            
            
            
            // writer.writeAdd(data, filesToMerge[i]->getEntryLen(entryId) - 1, 0);
        }
        for (auto it = memkeyset.begin(); it != memkeyset.end(); it++) {
            if (*it != repKey) {
                resultOutString.append(SSTR(*it));
                resultOutString.push_back('\n');
            }

        }
        writer.writeData(resultOutString.c_str(), resultOutString.size(), repKey,0);
        resultOutString.clear();
        // writer.writeEnd(key, 0);
    }
    writer.close();
    for (size_t i = 0; i < fileCount; i++) {
        filesToMerge[i]->close();
        delete filesToMerge[i];
    }
    delete[] filesToMerge;

    return EXIT_SUCCESS;
}
