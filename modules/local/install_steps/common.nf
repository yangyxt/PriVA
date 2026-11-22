// Shared by the install steps that run as their own process (registry
// "module" in conf/install_steps.json). Such a step's task gets
//
//   id, lineage    the step and its lineage (workflows/install.nf)
//   plan           plan.json
//   receipts       the receipts of the steps it needs, as inputs/N/receipt.json
//   planned        its plan entry: settings, release, pinned, user_paths
//   names          the file name in staged/ of each staged key
//   staged         the files it reads, staged into staged/: the keys of the
//                  steps it needs that its registry entry lists in "stage"
//                  (with their index files), and the files you supply for it
//
// and writes only into its work folder: its results into data/, then
// receipt.json from `priva install-task record`. storeDir moves both into
// resource_root/ID/LINEAGE when the task succeeds, and a later run finds them
// there and skips the task. These helpers write the record and adopt lines of
// the install runner (install-task), which no step's code digest counts; a
// step's own commands live in its module, whose text its digest counts.

// One shell word: the value in single quotes, each single quote in it written
// as '"'"'.
def installWord(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

// A process input that holds one file or several, as a list.
def installList(value) {
    value == null ? [] : (value instanceof List ? value : [value])
}

// The step's folder, where storeDir keeps its results.
def installFolder(String id, String lineage) {
    "${params.resource_root}/${id}/${lineage}"
}

// A setting of the step from the plan; a step stops when one it needs is unset.
def installSetting(planned, String key) {
    def value = planned.settings ? planned.settings[key] : null
    if (value == null || value.toString() == '' || value.toString() == 'null')
        throw new IllegalArgumentException("Set ${key} in the install configuration")
    value.toString()
}

// A setting of the step from the plan, or the fallback when it is unset.
def installValue(planned, String key, String fallback) {
    def value = planned.settings ? planned.settings[key] : null
    (value == null || value.toString() in ['', 'null']) ? fallback : value.toString()
}

// The commit a pinned git source was built from ("initial" before any
// refresh, or a stub run: then the source's current commit is taken).
def installCommit(planned) {
    def release = planned.release ? planned.release.toString() : ''
    if (!release.startsWith('[')) return ''
    def parts = new groovy.json.JsonSlurper().parseText(release)
    def commit = parts ? parts[0].release.toString() : ''
    commit ==~ /[0-9a-f]{40}/ ? commit : ''
}

// The first lines of a step's script: the priva program and the threads
// priva may use to hash the results.
def installStart(task) {
    """
    PRIVA_BIN=${installWord(params.priva_binary)}
    export INSTALL_THREADS=${task.cpus}
    mkdir -p data
    """.stripIndent()
}

// The last line of a step's script: priva checks the output keys (KEY: a path
// under data/, or a staged file you supply) and writes receipt.json.
def installRecord(String id, String lineage, receipts, Map keys) {
    def inputs = installList(receipts).collect { r -> "--input ${installWord(r)}" }.join(' ')
    def given = keys.collect { k, v -> "--key ${installWord(k + '=' + (v ?: ''))}" }.join(' ')
    "\"\$PRIVA_BIN\" install-task record --source-root ${installWord(projectDir)} --registry ${installWord("${projectDir}/conf/install_steps.json")} --plan plan.json --step ${id} --lineage ${lineage} --folder ${installWord(installFolder(id, lineage))} --receipt receipt.json ${inputs} ${given}"
}

// A step adopted from an existing installation (lineage adopted-...): only
// its receipt, naming that installation's files.
def installAdopt(String id, String lineage) {
    """
    mkdir -p data
    ${installWord(params.priva_binary)} install-task adopt --source-root ${installWord(projectDir)} --plan plan.json --step ${id} --lineage ${lineage} --receipt receipt.json
    """.stripIndent()
}

// The stub of every such step: a stub file for each output key.
def installStub(String id, String lineage, receipts, task) {
    if (lineage.startsWith('adopted-')) return installAdopt(id, lineage)
    installStart(task) + installRecord(id, lineage, receipts, [:]) + ' --stub\n'
}
