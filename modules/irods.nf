process get_file_from_irods { 
    tag "$sample_id"

    /**
    * Retrieves a file from iRODS with the given path
    */
    cache 'lenient'

    input:
        // declared as val to prevent stagging of an inexistent file
        tuple val(sample_id), val(irods_path)

    output:
        tuple val(sample_id), path(name)

    script:
        // if such a file exists then don't get it from iRODs
        // this assume that there is no collision between iRODs and file paths.
        // if we adopt a unabiguous naming scheme (e.g. add prefix irods:// or file://
        // we can get rid of such rare issue.
        name = file(irods_path).getName()
        
        """
        set -euxo pipefail
        ### for some reason without this line the HOME is set to /root despite singularity running
        ### using the invoker ordinary user:
        eval export HOME=~`whoami`
        # if target is a symlink we copy it so that we don't add another reference step
        if [[ -L ${irods_path} ]]; then
           // readlink is necessary to deal with possible relative symlinks:
            ln -s `readlink -f ${irods_path}` ${name} || exit 101
        elif [[ -f ${irods_path} ]]; then
            ln -s ${irods_path} ${name} || exit 102
        else 
            iget -K -f ${irods_path} -N ${task.cpus} ${name} || exit 103
            file_md5=`md5sum ${name} | awk '{print \$1}'` || exit 104
            irods_md5=`ichksum ${irods_path} | awk  '{print \$2}' | head -n 1` || exit 105
            if [[ "\$file_md5" = "\$irods_md5" ]]; then
                echo "MD5 matches"
            else 
                echo "MD5 mismatch: local md5 \$file_md5 != iRODs md5 \$irods_md5"
                exit 106
            fi
            # paranohia test as it has been shown that irods recovered files had some changes
            samtools quickcheck -u ${name} || exit 106
            samtools view -bS ${name} > /dev/null || exit 107
        fi
        """ 
}