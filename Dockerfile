FROM julia
COPY ./Project.toml /app/Project.toml
COPY ./Manifest.toml /app/Manifest.toml
COPY ./src /app/src
COPY ./test /app/test

WORKDIR /app

RUN apt-get update && apt-get install -y \
    cmake \
    xorg-dev \
    mesa-utils \
    xvfb \
    libgl1 \
    freeglut3-dev

RUN julia -e 'using Pkg; Pkg.add("Revise");'
RUN julia -e 'using Pkg; Pkg.activate("./"); Pkg.instantiate();'

CMD ["julia", "--project", "/app", "-e", "using CylindersBasedCameraResectioning; main()"]
