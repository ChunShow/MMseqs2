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
#include <vector>

struct CountTableEntry {
    unsigned int proteinId;
    unsigned int entryCount, clusterSizeWeightedCount, seqLength;
    CountTableEntry() : proteinId(UINT_MAX), entryCount(0), clusterSizeWeightedCount(0), seqLength(0) {}
    CountTableEntry(unsigned int proteinId, unsigned int entryCount, unsigned int clusterSizeWeightedCount, unsigned int seqLength) : proteinId(proteinId), entryCount(entryCount), clusterSizeWeightedCount(clusterSizeWeightedCount), seqLength(seqLength) {}
};

struct DBEntry {
    unsigned int repId;
    unsigned int proteinId;
    DBEntry() : repId(UINT_MAX) {}
    DBEntry(unsigned int repId, unsigned int proteinId) : repId(repId), proteinId(proteinId) {}
    
    static bool compareByRepIdNProteinId(const DBEntry &a, const DBEntry &b)
    {
        if (a.repId < b.repId) {
            return true;
        }
        if (a.repId > b.repId) {
            return false;
        }

        bool aIsRep = (a.repId == a.proteinId); 
        bool bIsRep = (b.repId == b.proteinId);
        
        // move Rep to the front
        if (aIsRep && !bIsRep) {
            return true;    
        }
        if (!aIsRep && bIsRep) {
            return false; 
        }

        return (a.proteinId < b.proteinId);
    }
    
    static bool compareByCountTable(const DBEntry &a, const DBEntry &b, 
                                const std::vector<CountTableEntry> &countTable)
    {

        const bool aIsPrevRep = (a.repId == a.proteinId);
        const bool bIsPrevRep = (b.repId == b.proteinId);
        
        // If a or b is prevRep, it should be the last in group
        if (aIsPrevRep && !bIsPrevRep) {
            return false;  
        }
        if (!aIsPrevRep && bIsPrevRep) {
            return true;   
        }

        const CountTableEntry &ca = countTable[a.proteinId];
        const CountTableEntry &cb = countTable[b.proteinId];
        // Priority: 
        // 1) entryCount 
        if (ca.entryCount > cb.entryCount) return true;
        if (ca.entryCount < cb.entryCount) return false;
        
        // 2) sequenceLength 
        if (ca.seqLength > cb.seqLength) return true;
        if (ca.seqLength < cb.seqLength) return false;

        // // 1) clusterSizeWeightedCount 
        // if (ca.clusterSizeWeightedCount > cb.clusterSizeWeightedCount) return true;
        // if (ca.clusterSizeWeightedCount < cb.clusterSizeWeightedCount) return false;

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


    std::vector<CountTableEntry> countTable;
    countTable.reserve(qdbr.getSize());
    std::cout << "qdbr size: " << qdbr.getSize() << std::endl;
    for (size_t i=0; i < qdbr.getSize(); ++i) {
        unsigned int dbKey = qdbr.getDbKey(i);
        size_t queryId = qdbr.getId(dbKey);
        
        unsigned int seqLength = qdbr.getSeqLen(queryId);
        countTable.push_back(CountTableEntry(queryId, 0, 0, 0));
    }

    if (countTable.size() != qdbr.getSize()) {
        std::cout << "Error: " << countTable.size() << " " << qdbr.getSize() << std::endl;
        EXIT(EXIT_FAILURE);
    }
    std::vector<DBEntry> DBEntries;
    DBEntries.reserve(resdbr.getSize());

    #pragma omp parallel 
    {
        unsigned int thread_idx = 0;
    #ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
    #endif    
        char buffer[1024 + 32768 * 4];
        std::vector<DBEntry> localDBEntries;
        #pragma omp for schedule(dynamic, 1)
        for (size_t id=0; id < resdbr.getSize(); ++id) {
            char* data = resdbr.getData(id, thread_idx);
            std::vector<unsigned int> memberIds;
            //1. Fill countTable
            unsigned int repId = UINT_MAX;
            if (*data != '\0') {
                Util::parseKey(data, buffer);
                unsigned int referenceKey = (unsigned int) strtoul(buffer, NULL, 10);
                repId= qdbr.getId(referenceKey);
            }
            while (*data != '\0') {
                Util::parseKey(data, buffer);
                const unsigned int key = (unsigned int) strtoul(buffer, NULL, 10);
                unsigned int proteinId= qdbr.getId(key);
                memberIds.push_back(proteinId);
                localDBEntries.push_back(DBEntry(repId, proteinId));
                data = Util::skipLine(data);
            }
            
            if (repId != memberIds[0]) {
                std::cout << "Error: " << repId << " " << memberIds[0] << std::endl;
                EXIT(EXIT_FAILURE);
            }

            unsigned int clustSize = memberIds.size();
            if (memberIds.size() > 1) { // Try to ignore singletons
                for (size_t i=0; i < memberIds.size(); ++i) {
                    __sync_fetch_and_add(&countTable[memberIds[i]].entryCount, 1);
                    __sync_fetch_and_add(&countTable[memberIds[i]].clusterSizeWeightedCount, clustSize);
                }
            }
        }

        #pragma omp critical
        {
            DBEntries.insert(DBEntries.end(),
                            std::make_move_iterator(localDBEntries.begin()),
                            std::make_move_iterator(localDBEntries.end()));
        }
    }

    // load DBEntries and sort its members by countTable
    unsigned int prevRepId = DBEntries[0].repId;
    size_t prevGroupStart = 0;
    for (size_t elementIdx = 0; elementIdx <= DBEntries.size(); ++elementIdx) { 
        unsigned int currRepId = (elementIdx == DBEntries.size()) ? UINT_MAX : DBEntries[elementIdx].repId;
        
        if (prevRepId != currRepId) {
            // if singletons, do nothing
            if (elementIdx == prevGroupStart + 1) {
                prevRepId = currRepId;
                prevGroupStart = elementIdx;
                continue;
            }
            // else sort the group by countTable
            SORT_SERIAL(DBEntries.begin() + prevGroupStart, DBEntries.begin() + elementIdx,
                    [&](const DBEntry &lentry, const DBEntry &rentry) {
                        return DBEntry::compareByCountTable(lentry, rentry, countTable);
                    });

            unsigned int newRepId = DBEntries[prevGroupStart].proteinId;
            // reassign(change) repId
            for (size_t i = prevGroupStart; i < elementIdx; ++i) {
                DBEntries[i].repId = newRepId;
            }
            prevRepId = currRepId;
            prevGroupStart = elementIdx;
        }
    }
    // sort DBEntries by new RepKey
    // 1) sort by repId 2) In each group, sort by proteinId. Rep should be the first in the group for writer
    SORT_SERIAL(DBEntries.begin(), DBEntries.end(), DBEntry::compareByRepIdNProteinId);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_PREFILTER_RES);
    dbw.open();
    
    prevRepId = DBEntries[0].repId;
    prevGroupStart = 0;
    if (prevRepId != 0 ) {
        //write singleton result
        std::string resultOutString;
        resultOutString.reserve(300);
        for (size_t id = 0; id < prevRepId; ++id) {
            resultOutString.append(SSTR(qdbr.getDbKey(id)));
            resultOutString.push_back('\n');
            dbw.writeData(resultOutString.c_str(), resultOutString.length(), qdbr.getDbKey(id), 0);
            resultOutString.clear();
        }
    }

    for (size_t elementIdx = 0; elementIdx <= DBEntries.size(); ++elementIdx) {
        unsigned int currRepId = (elementIdx == DBEntries.size()) ? qdbr.getSize() : DBEntries[elementIdx].repId; // Need to check 
        if (currRepId != prevRepId) {
            std::string resultOutString;
            resultOutString.reserve(1024);

            // write the current group
            unsigned int prevProtId = UINT_MAX;
            for (size_t i = prevGroupStart; i < elementIdx; ++i) { //instead of using set
                if (prevProtId != DBEntries[i].proteinId) {
                    resultOutString.append(SSTR(qdbr.getDbKey(DBEntries[i].proteinId)));
                    resultOutString.push_back('\n');
                    prevProtId = DBEntries[i].proteinId;
                }
            }
            dbw.writeData(resultOutString.c_str(), resultOutString.length(), qdbr.getDbKey(prevRepId), 0);
            resultOutString.clear();
            //remaining singletons : prevRepId + 1 ~ currRepId - 1
            if (currRepId - prevRepId > 1) { // discontinuous
                for (size_t rep = prevRepId + 1; rep < currRepId; ++rep) {
                    resultOutString.append(SSTR(qdbr.getDbKey(rep)));
                    resultOutString.push_back('\n');
                    dbw.writeData(resultOutString.c_str(), resultOutString.length(), qdbr.getDbKey(rep), 0);
                    resultOutString.clear();
                }
            }
            prevRepId = currRepId;
            prevGroupStart = elementIdx;
        }
    }

    dbw.close(true);
    qdbr.close();
    resdbr.close();
    return EXIT_SUCCESS;
}