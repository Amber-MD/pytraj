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
                sh "conda install -y --file conda-requirements.txt"
                sh "pip install -r pip-requirements.txt"
                sh "python setup.py install --user"
                sh "python run_tests.py"
            }
        }
    }
}
