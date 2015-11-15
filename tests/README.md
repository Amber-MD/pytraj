- Run test

    - nosetests --with-coverage --cover-package pytraj -vs .
    - nosetests -vs --processes 4 --process-timeout 200 .
    - nosetests -vs --processes 4 --process-timeout 200 --with-coverage --cover-package pytraj .

- Format

    - yapf --style pep8 ./tests/*py tests/*/*py pytraj/*py pytraj/*/*py
    - autopep8 -a my_file.py -i
    - eradicate -i my_file.py 
