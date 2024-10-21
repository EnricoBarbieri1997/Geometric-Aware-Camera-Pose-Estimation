FROM julia
COPY ./Project.toml /app/Project.toml
COPY ./Manifest.toml /app/Manifest.toml
COPY ./src /app/src
WORKDIR /app
RUN julia -e 'using Pkg; Pkg.add("Revise");'
RUN julia -e 'using Pkg; Pkg.activate("./"); Pkg.instantiate();'

CMD ["julia", "-e", "using CylindersBasedCameraResectioning; main()"]
