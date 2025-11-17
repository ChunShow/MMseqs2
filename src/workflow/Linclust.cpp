#include "Parameters.h"
#include "Util.h"
#include "DBWriter.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"

#include "linclust.sh.h"

#include <cassert>

void setLinclustWorkflowDefaults(Parameters *p) {
    p->spacedKmer = false;
    p->covThr = 0.8;
    p->maskMode = 0;
    p->evalThr = 0.001;
    p->seqIdThr = 0.9;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID; // use alignment mode 3 in linclust2
    p->skipHamming = true;
}

int linclust(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    setLinclustWorkflowDefaults(&par);
    par.PARAM_ADD_BACKTRACE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ALT_ALIGNMENT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_ZDROP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_RESCORE_MODE.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_MAX_ACCEPT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.overrideParameterDescription(par.PARAM_S, "Sensitivity will be automatically determined but can be adjusted", NULL, par.PARAM_S.category | MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_INCLUDE_ONLY_EXTENDABLE.addCategory(MMseqsParameter::COMMAND_EXPERT);

    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.linclustworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());

    // # 1. Finding exact $k$-mer matches.
    bool kmerSizeWasSet = false;
    bool alphabetSizeWasSet = false;
    bool clusterModeSet = false;
    for (size_t i = 0; i < par.linclustworkflow.size(); i++) {
        if (par.linclustworkflow[i]->uniqid == par.PARAM_K.uniqid && par.linclustworkflow[i]->wasSet) {
            kmerSizeWasSet = true;
        }
        if (par.linclustworkflow[i]->uniqid == par.PARAM_ALPH_SIZE.uniqid && par.linclustworkflow[i]->wasSet) {
            alphabetSizeWasSet = true;
        }
        if (par.linclustworkflow[i]->uniqid == par.PARAM_CLUSTER_MODE.uniqid && par.linclustworkflow[i]->wasSet) {
            clusterModeSet = true;
        }
    }

    const bool nonSymetric = (par.covMode == Parameters::COV_MODE_TARGET || par.covMode == Parameters::COV_MODE_QUERY);
    if (clusterModeSet == false){
        if (nonSymetric) {
            par.clusteringMode = Parameters::GREEDY_MEM;
        } else {
            par.clusteringMode = Parameters::SET_COVER;
        }
        std::string cluMode = (par.clusteringMode==Parameters::GREEDY_MEM) ? "GREEDY MEM" : "SET COVER";
        Debug(Debug::INFO) << "Set cluster mode " << cluMode << ".\n";
    }

    if (kmerSizeWasSet == false) {
        par.kmerSize = Parameters::CLUST_LINEAR_DEFAULT_K;
    }
    if (alphabetSizeWasSet == false) {
        par.alphabetSize = MultiParam<NuclAA<int>>(NuclAA<int>(Parameters::CLUST_LINEAR_DEFAULT_ALPH_SIZE, 5));
    }

    const int dbType = FileUtil::parseDbType(par.db1.c_str());
    const bool isUngappedMode = par.alignmentMode == Parameters::ALIGNMENT_MODE_UNGAPPED;
    if (isUngappedMode && Parameters::isEqualDbtype(dbType, Parameters::DBTYPE_HMM_PROFILE)) {
        par.printUsageMessage(command, MMseqsParameter::COMMAND_ALIGN|MMseqsParameter::COMMAND_PREFILTER);
        Debug(Debug::ERROR) << "Cannot use ungapped alignment mode with profile databases.\n";
        EXIT(EXIT_FAILURE);
    }

    
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());
    cmd.addVariable("ALIGNBLOCK_PAR", par.createParameterString(par.alignblock).c_str());
    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());
    std::string program = tmpDir + "/linclust.sh";
    FileUtil::writeFile(program, linclust_sh, linclust_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Unreachable
    assert(false);
    return 0;
}