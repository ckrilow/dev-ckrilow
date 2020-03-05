# Pipeline: example_module_1
# Author: Forename Surname
# Summary: Does something

workflow example_module_1 {
    # Pipeline version
    String pipeline_ver = "0.1.0"

    # Parameters ###############################################################
    # Initialize required parameters
    File my_file1                           # input file stat file from
                                            # blah pipeline.
                                            # Note: should be the full path

    # Initialize parameters with defaults
    String scripts_directory = "/bin"       # directory containing all of the
                                            # scripts for this workflow
    String output_file_tag = "my_workflow_results"
                                            # tag for output file.
    String filesystem = "local"             # Valid options: [local, google].
                                            # Specifies method to (a) see if a
                                            # file exists and (b) to copy the
                                            # results to the root dir
    Boolean force_run_all_tasks = false     # If true, all steps of the workflow
                                            # are run. If false, the workflow
                                            # checks to see if the COJO files
                                            # are present in skips generating
                                            # them if so.

    # Resources
    Int threads = "16"                      # threads
    String temp_dir = "/tmp"                # temp dir where temp files will be
                                            # held and then deleted
    ############################################################################

    # Begin pipeline ###########################################################
    # get the path to the input in order to copy key files there
    # the input_file_dir path variable will end in a backslash
    String input_file_dir = sub(
        my_file1,
        basename(my_file1),
        ""
    )
    String out_dir_final_results = input_file_dir + "example_module_1"
    String out_results = out_dir_final_results + "/" +
        output_file_tag + "-example_module_1.tsv.gz"


    # if force_run_all_tasks == True then always call runAnalysis.
    # if not force_run_all_tasks, then look at the input directory to see
    # if there are output files
    if (!force_run_all_tasks) {
        Array[String] files_to_check = [out_results]
        call doFilesExist {
            input:
                path_to_files=files_to_check,
                filesystem=filesystem
        }
    }
    Boolean run_analysis = if (force_run_all_tasks == true)
        then true
        else (doFilesExist.out == false)


    if (run_analysis == true) {
        call runAnalysis {
            input:
                my_file=my_file1,
                output_file_tag=output_file_tag,
                scripts_directory=scripts_directory,
                verbose=true,
                threads=threads
        }

        # copy the final file back to the root directory
        Array[File] final_results_files = [runAnalysis.results]
        call copyFiles as copyFilesFinalResults {
            input:
                files_to_copy=final_results_files,
                destination=out_dir_final_results,
                filesystem=filesystem
        }
    }
    # set files for downstream tasks based on if we re-made them.
    File? results = if (run_analysis == true)
        then runAnalysis.results
        else out_results

    output {
        File? results = results
    }
    meta {
        author: "Forename Surname"
    }
} # end workflow::example_module_1


task printString {
    String input_string
    command {
        echo -e "input_string:\t${input_string}"
    }
    output {
    	File response = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
        cpu: "1"
        memory: "500 MB"
        disks: "local-disk 1 SSD" # minimum size: local-disk 1 SSD
    }
}


task doFilesExist {
    Array[String] path_to_files
    String filesystem

    # for a single file, the below code works where:
    # 0 means exist and 1 means does not exist
    # gsutil -q stat gs://main-bucket/sub-directory-bucket/*; echo $?
    # (ls ${path_to_file} >> /dev/null 2>&1 && echo 0) || echo 1
    command <<<
        python <<CODE
        import os
        from google.cloud import storage
        from urllib.parse import urlparse
        def gcs_does_file_exist(filepath):
            p = urlparse(filepath)
            client = storage.Client()
            bucket = client.get_bucket(p.netloc)
            blob = bucket.blob(p.path.strip("/")) # can't start with backslah
            return blob.exists()
        files = "${sep=',' path_to_files}".split(",")
        all_files_exist = True
        for f in files:
            if "${filesystem}" == "local":
                if not os.path.isfile(f):
                    all_files_exist = False
                    break
            elif "${filesystem}" == "google":
                if not gcs_does_file_exist(f):
                    all_files_exist = False
                    break
            else:
                Exception(
                    "ERROR: filesystem=${filesystem} valid parameters are",
                    " [local, google]."
                )
        print(str(all_files_exist).lower())
        CODE
    >>>
    output {
        Boolean out = read_boolean(stdout())
    }
    runtime {
        docker: "python:latest"
        cpu: "1"
        memory: "1024 MB"
        disks: "local-disk 1 SSD" # minimum size: local-disk 1 SSD
    }
}


task copyFiles {
    Array[File] files_to_copy
    String destination
    String filesystem

    # Copy files differently based on if the file location is local or on
    # google cloud.
    command {
        if [ "${filesystem}" == "local" ]
        then
            mkdir -p ${destination}
            cp -L -R -u ${sep=' ' files_to_copy} ${destination}
        elif [ "${filesystem}" == "google" ]
        then
            gsutil cp -I ${destination} < ${write_lines(files_to_copy)}
        else
            (>&2 echo "ERROR: filesystem=${filesystem} valid parameters are [local, google].")
            exit 1
        fi
    }
    runtime {
        docker: "google/cloud-sdk"
        cpu: "1"
        memory: "500 MB"
        disks: "local-disk 1 SSD" # minimum size: local-disk 1 SSD
    }
}


task cleanUpFiles {
  Array[File] files_to_delete
  String filesystem

  command {
      if [ "${filesystem}" == "local" ]
      then
          rm -r -f ${sep=' ' files_to_delete}
      elif [ "${filesystem}" == "google" ]
      then
          gsutil rm -I < ${write_lines(files_to_delete)}
      else
          (>&2 echo "ERROR: filesystem=${filesystem} valid parameters are [local, google].")
          exit 1
      fi
  }
  runtime {
      docker: "google/cloud-sdk"
      cpu: "1"
      memory: "1024 MB"
      disks: "local-disk 1 SSD" # minimum size: local-disk 1 SSD
  }
}


task runAnalysis {
    # Runs analysis

    # required parameters:
    File my_file                            # input_file

    # parameters with defaults
    String? file_id_bait                    # descriptor here
    String output_file_tag = "results"
                                            # tag for output file.
    String scripts_directory = "/bin"       # directory containing all of the
                                            # scripts for this workflow
    Boolean verbose = false                 # run in verbose mode
    Int threads = "2"                       # threads

    # private variables
    File script__run_script = scripts_directory + "/my_script.R"
    Int runtime_cpu = ceil(threads / 2)

    command {
        Rscript ${script__run_script} \
            ${"--file " + my_file} \
            ${"--threads " + threads} \
            ${if verbose then "--verbose " else ""}
    }
    output {
        File results = "${output_file_tag}-results.tsv.gz"
    }
    runtime {
        docker: "us.gcr.io/foo:latest"
        cpu: runtime_cpu
        memory: "25240 MB"
        disks: "local-disk 18 SSD"
        failOnStderr: false
    }
}
