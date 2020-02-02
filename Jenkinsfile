#!/usr/bin/env groovy

pipeline {
    agent none

    parameters {
        string defaultValue: "master", description: "The libcpptraj docker image tag to use",
               name: "LIBCPPTRAJ_IMAGE_TAG", trim: true
        string defaultValue: "master", description: "The pytraj commit/branch to test",
               name: 'BRANCH_TO_BUILD', trim: true
    }

    post {
        failure {
            emailext attachLog: true, compressLog: true,
                     subject: "pytraj tests failed",
                     body: "The pytraj tests failed when running against the ambermd/libcpptraj:${LIBCPPTRAJ_IMAGE_TAG} docker image",
                     recipientProviders: [culprits(), brokenTestsSuspects(), developers()]
        }
    }

    stages {
        stage("Build and test pytraj") {
            agent {
                docker{ image "ambermd/libcpptraj:${LIBCPPTRAJ_IMAGE_TAG}" }
            }

            steps {
                sh "pip install -r pip-requirements.txt"
                sh "python setup.py install --user"
                // 20 minute timeout, then kill the job. If it doesn't hang, it works in <10min
                timeout(20) {
                    sh "cd tests && pytest -vs --ignore=test_parallel_pmap --ignore=test_run_mpi.py --ignore=test_energy/test_pmap_sander.py --ignore=test_parallel_mpi --ignore=test_actionlist.py"
                }
            }
        }
    }
}
