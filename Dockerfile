FROM julia:1.11.5-bookworm
COPY ./Project.toml /app/Project.toml
COPY ./Manifest.toml /app/Manifest.toml
COPY ./src /app/src
COPY ./test /app/test
COPY ./scripts /app/scripts

WORKDIR /app

RUN apt-get update && apt-get install -y \
    cmake \
    xorg-dev \
    mesa-utils \
    xvfb \
    libgl1 \
    freeglut3-dev \
    inotify-tools

RUN mkdir ~/.julia
RUN mkdir ~/.julia/config
RUN printf "atreplinit() do repl\ntry\n@eval using Revise\n@async Revise.wait_steal_repl_backend()\ncatch\nend\nend\n" >> ~/.julia/config/startup.jl 

RUN julia --project=. -e 'using Pkg; Pkg.add("Revise");'
RUN julia --project=. -e 'using Pkg; Pkg.update(); Pkg.resolve(); Pkg.instantiate();'

CMD ["julia", "--project", "/app", "-e", "using CylindersBasedCameraResectioning; main()"]
