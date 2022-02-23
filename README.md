# boltz-t-cpp

## Установка

### Docker
1. Установите Docker https://www.docker.com/
2. В командной строке выполните docker ```pull intel/oneapi-hpckit```
эта команда скачивает Docker image с установленными Intel библиотеками и компиляторами.
3. Скачайте код с репозитория. В командной строке выполните
```docker run -t --name boltz-cpp -v PATH_TO_REPO:/repo intel/oneapi-hpckit:latest```
эта команда создаёт Docker container на основе image intel/oneapi-hpckit:latest
"-t" открывает командную строку внутри контейнера (может глючить), "--name" - задаёт имя контейнера, "-v" создаёт внутри контейнера папку ```/repo```, которая ссылается на папку ```PATH_TO_REPO``` в файловой системе.
4. В командной строке контейнера (она доступна также в GUI Docker) перейдите в папку /repo/src и выполните ```make```

### WSL (Windows)

### Ручная установка



