process cleanup {
    executor = 'local'
    container false
    
    cpus params.CPU0
    memory params.MEM0 + "g"
    clusterOptions = "-clear -cwd -l vf=${memory},num_proc=${cpus} -binding linear:${cpus} " + (params.project.equalsIgnoreCase("none")? "" : "-P " + params.project) + " -q ${params.queue} ${params.extraCluOpt}"

    when:
    !params.keepFiles

    input:
    path(files)

    script:
    println files
    // Use 'true' to always succeed even if the file doesn't exist or can't be deleted
    """
    rm -f `realpath $files` || true
    """
}
