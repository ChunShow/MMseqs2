#include "Util.h"
#include "Parameters.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "MemoryMapped.h"
#include "Alignment.h"
#include "itoa.h"
#include "Timer.h"
#include <algorithm>
#include <unordered_set>
#include <vector>

struct CountEntry {
    unsigned int proteinKey;
    unsigned int entryCount, clusterSizeWeightedCount, repWeightedCount;
    CountEntry() : proteinKey(UINT_MAX), entryCount(0), clusterSizeWeightedCount(0), repWeightedCount(0) {}
    CountEntry(unsigned int proteinKey, unsigned int entryCount, unsigned int clusterSizeWeightedCount, unsigned int repWeightedCount) : proteinKey(proteinKey), entryCount(entryCount), clusterSizeWeightedCount(clusterSizeWeightedCount), repWeightedCount(repWeightedCount) {}
};

struct prefEntry {
    unsigned int repKey;
    unsigned int proteinKey;
    prefEntry() : repKey(UINT_MAX) {}
    prefEntry(unsigned int repKey, unsigned int proteinKey) : repKey(repKey), proteinKey(proteinKey) {}
    
    static bool compareByRepKeyNProteinKey(const prefEntry &a, const prefEntry &b)
    {
        if (a.repKey < b.repKey) {
            return true;
        }
        if (a.repKey > b.repKey) {
            return false;
        }

        bool a_special = (a.repKey == a.proteinKey); 
        bool b_special = (b.repKey == b.proteinKey);

        if (a_special && !b_special) {
            return true;    
        }
        if (!a_special && b_special) {
            return false; 
        }

        return (a.proteinKey < b.proteinKey);
    }
    static bool compareByCountTable(const prefEntry &a, const prefEntry &b, 
                                const std::vector<CountEntry> &countTable)
    {

        // const bool a_is_special = (a.repKey == a.proteinKey);
        // const bool b_is_special = (b.repKey == b.proteinKey);

        // if (a_is_special && !b_is_special) {
        //     return false;  
        // }
        // if (!a_is_special && b_is_special) {
        //     return true;   
        // }

        const CountEntry &ca = countTable[a.proteinKey];
        const CountEntry &cb = countTable[b.proteinKey];
        
      
        // 2) entryCount 
        if (ca.entryCount > cb.entryCount) return true;
        if (ca.entryCount < cb.entryCount) return false;

        // 1) clusterSizeWeightedCount 
        if (ca.clusterSizeWeightedCount > cb.clusterSizeWeightedCount) return true;
        if (ca.clusterSizeWeightedCount < cb.clusterSizeWeightedCount) return false;


        // // 3) repWeightedCount 
        // if (ca.repWeightedCount > cb.repWeightedCount) return true;
        // if (ca.repWeightedCount < cb.repWeightedCount) return false;
       
        return false;
    }
};


int resortprefilter(int argc, const char **argv, const Command &command){
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_LOOKUP);
    qdbr.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> resdbr(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    resdbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);


    std::vector<CountEntry> countTable; // need to optimize
    countTable.reserve(qdbr.getSize());

    for (size_t i=0; i < qdbr.getSize(); ++i) {
        unsigned int dbKey = qdbr.getDbKey(i);
        countTable.push_back(CountEntry(dbKey, 0, 0, 0));
        // std::cout << i << " " << dbKey << std::endl;
        // if (i != dbKey) {
        //     std::cout << "Error: " << i << " " << dbKey << std::endl;
        // }
    }

    if (countTable.size() != qdbr.getSize()) {
        std::cout << "Error: " << countTable.size() << " " << qdbr.getSize() << std::endl;
        EXIT(EXIT_FAILURE);
    }

    std::vector<prefEntry> prefEntries;
    prefEntries.reserve(resdbr.getSize());

    #pragma omp parallel 
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
    #endif    
        char buffer[1024 + 32768 * 4];
        std::vector<prefEntry> localPrefEntries;
        #pragma omp for schedule(dynamic, 1)
        for (size_t id=0; id < resdbr.getSize(); ++id) {
            char* data = resdbr.getData(id, thread_idx);
            std::vector<unsigned int> memberKeys;
            //1. Fill countTable
            unsigned int referenceKey = UINT_MAX;
            if (*data != '\0') {
                Util::parseKey(data, buffer);
                referenceKey = (unsigned int) strtoul(buffer, NULL, 10);
            }
            while (*data != '\0') {
                Util::parseKey(data, buffer);
                const unsigned int key = (unsigned int) strtoul(buffer, NULL, 10);
                memberKeys.push_back(key);
                localPrefEntries.push_back(prefEntry(referenceKey, key));
                data = Util::skipLine(data);
            }
            
            if (referenceKey != memberKeys[0]) {
                std::cout << "Error: " << referenceKey << " " << memberKeys[0] << std::endl;
                EXIT(EXIT_FAILURE);
            }

            unsigned int clustSize = memberKeys.size();
            if (memberKeys.size() > 1) { // Try to ignore singletons
                for (size_t i=0; i < memberKeys.size(); ++i) {
                    __sync_fetch_and_add(&countTable[memberKeys[i]].entryCount, 1); // update EntryCount
                    __sync_fetch_and_add(&countTable[memberKeys[i]].clusterSizeWeightedCount, clustSize);
                    __sync_fetch_and_add(&countTable[memberKeys[i]].repWeightedCount, 1);
                }
                __sync_fetch_and_add(&countTable[referenceKey].repWeightedCount, clustSize-1);
            }
        }

        #pragma omp critical
        {
            prefEntries.insert(prefEntries.end(),
                            std::make_move_iterator(localPrefEntries.begin()),
                            std::make_move_iterator(localPrefEntries.end()));
        }
    }

    // load prefEntries and sort its members by countTable
    unsigned int prevRepKey = prefEntries[0].repKey;
    size_t prevGroupStart = 0;
    for (size_t elementIdx = 0; elementIdx <= prefEntries.size(); ++elementIdx) { //등호?
        unsigned int currRepKey = (elementIdx == prefEntries.size()) ? UINT_MAX : prefEntries[elementIdx].repKey;
        
        if (prevRepKey != currRepKey) {
            // if singletons, do nothing
            if (elementIdx == prevGroupStart + 1) {
                prevRepKey = currRepKey;
                prevGroupStart = elementIdx;
                continue;
            }
            // std::cout << "group: " << prevRepKey << std::endl;
            // else sort the group and reset repKey
            //Debug before sort
            SORT_SERIAL(prefEntries.begin() + prevGroupStart, prefEntries.begin() + elementIdx,
                    [&](const prefEntry &lhs, const prefEntry &rhs) {
                        return prefEntry::compareByCountTable(lhs, rhs, countTable);
                    });

            unsigned int newRepKey = prefEntries[prevGroupStart].proteinKey;
            //Resetup repKey
            for (size_t i = prevGroupStart; i < elementIdx; ++i) {
                prefEntries[i].repKey = newRepKey;
            }
            prevRepKey = currRepKey;
            prevGroupStart = elementIdx;
        }

    }
    SORT_SERIAL(prefEntries.begin(), prefEntries.end(), prefEntry::compareByRepKeyNProteinKey);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    dbw.open();
    
    prevRepKey = prefEntries[0].repKey;
    prevGroupStart = 0;
    if (prevRepKey != 0 ) {
        //write singleton result
        std::string resultOutString;
        resultOutString.reserve(300);
        for (size_t key = 0; key < prevRepKey; ++key) {
            resultOutString.append(SSTR(key));
            resultOutString.push_back('\n');
            dbw.writeData(resultOutString.c_str(), resultOutString.length(), key, 0);
            resultOutString.clear();
        }
    }

    for (size_t elementIdx = 0; elementIdx <= prefEntries.size(); ++elementIdx) {
        unsigned int currRepKey = (elementIdx == prefEntries.size()) ? qdbr.getSize() : prefEntries[elementIdx].repKey; // Need to check 
        if (currRepKey != prevRepKey) {
            std::string resultOutString;
            resultOutString.reserve(1024);

            // write the current group
            unsigned int prevProtKey = UINT_MAX;
            for (size_t i = prevGroupStart; i < elementIdx; ++i) { //instead of using set
                if (prevProtKey != prefEntries[i].proteinKey) {
                    resultOutString.append(SSTR(prefEntries[i].proteinKey));
                    resultOutString.push_back('\n');
                    prevProtKey = prefEntries[i].proteinKey;
                }
            }
            dbw.writeData(resultOutString.c_str(), resultOutString.length(), prevRepKey, 0);
            resultOutString.clear();
            //remaining singletons : prevRepKey + 1 ~ currRepKey - 1
            if (currRepKey - prevRepKey > 1) { // discontinuous
                for (size_t rep = prevRepKey + 1; rep < currRepKey; ++rep) {
                    resultOutString.append(SSTR(rep));
                    resultOutString.push_back('\n');
                    dbw.writeData(resultOutString.c_str(), resultOutString.length(), rep, 0);
                    resultOutString.clear();
                }
            }
            prevRepKey = currRepKey;
            prevGroupStart = elementIdx;
        }
    }

    dbw.close(true);
    qdbr.close();
    resdbr.close();
    return EXIT_SUCCESS;
}