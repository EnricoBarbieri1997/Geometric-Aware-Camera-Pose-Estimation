# PowerShell script to run Julia in Docker with different noise levels

# Define the noise levels
$noiseLevels = @()
for ($n = 0.00; $n -le 0.04; $n += 0.0016) {
    $noiseLevels += [math]::Round($n, 6)
}

# Path to your project directory (one level up from this script)
$projectPath = Resolve-Path "$PSScriptRoot\.."

# Docker image name
$imageName = "julia-cylinders"

foreach ($noise in $noiseLevels) {
    Write-Host "Running with noise level: $noise"

    docker run --rm `
        -v "$($projectPath):/app" `
        -v "$($projectPath)/juliacache:/root/.julia" `
        -w /app `
        --gpus all `
        -e DISPLAY=":0" `
        -e WAYLAND_DISPLAY="wayland-0" `
        -e XDG_RUNTIME_DIR="/mnt/wslg/runtime-dir" `
        -e PULSE_SERVER="/mnt/wslg/PulseServer" `
        -v "\\wsl.localhost\Ubuntu-18.04\tmp\.X11-unix:/tmp/.X11-unix" `
        -v "\\wsl.localhost\Ubuntu-18.04\mnt\wslg:/mnt/wslg" `
        $imageName `
        julia --project=. scripts/run_noise.jl $noise
}
