executable = dataGeneration.sh
getenv = True
arguments = "$(CLUSTER) $(PROCESS)"
output = $OUTPUT_DIR/outputfile.$(CLUSTER).$(PROCESS)
error = $OUTPUT_DIR/errorfile.$(CLUSTER).$(PROCESS)
log = $OUTPUT_DIR/example.job.$(CLUSTER).$(PROCESS).log
+MaxRuntime = {{ max_runtime }}
queue
