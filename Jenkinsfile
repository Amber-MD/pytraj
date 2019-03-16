pipeline {
    agent none

    parameters {
        string defaultValue: "master", description: "The libcpptraj docker image tag to use",
               name: "LIBCPPTRAJ_IMAGE_TAG", trim: true
    }

    post {
        failure {
            emailext attachLog: true, compressLog: true,
                     subject: "pytraj tests failed",
                     body: "The pytraj tests failed when running against the ambermd/libcpptraj:${LIBCPPTRAJ_IMAGE_TAG}",
                     recipientProviders: [culprits(), brokenTestsSuspects(), developers()]
        }
    }

    stages {
        stage("Build and test pytraj") {
            agent {
                docker{ image "ambermd/libcpptraj:${LIBCPPTRAJ_IMAGE_TAG}" }
            }

            steps {
                sh "python setup.py install --user"
                sh "python run_tests.py --simple"
            }
        }
    }
}