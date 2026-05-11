executable = dataGeneration.sh
getenv = True
OUTPUT_DIR = $ENV(OUTPUT_DIR)
infile = $ENV(infile)
arguments = "$(OUTPUT_DIR) $(CLUSTER) $(PROCESS)"
output = $(OUTPUT_DIR)/jobOutput/outputfile.$(CLUSTER).$(PROCESS)
error = $(OUTPUT_DIR)/jobOutput/errorfile.$(CLUSTER).$(PROCESS)
log = $(OUTPUT_DIR)/jobOutput/job.$(CLUSTER).$(PROCESS).log
+MaxRuntime = {{ max_runtime }}
queue
