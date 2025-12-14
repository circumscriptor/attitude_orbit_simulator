# Approximate TNS-0 Simulation
./build/attitude_orbit_simulator \
    --output tns0_validation.csv \
    --mass 5000 \
    --width 0.2 --height 0.2 --length 0.2 \
    --magnet-remanence 1.0 \
    --magnet-diameter 0.01 --magnet-length 0.03 \
    --rod-volume 2.4e-7 \
    --rod-orientation 1,0,0 --rod-orientation 1,0,0 --rod-orientation 1,0,0 --rod-orientation 1,0,0 \
    --rod-orientation 0,1,0 --rod-orientation 0,1,0 --rod-orientation 0,1,0 --rod-orientation 0,1,0 \
    --hysteresis-ms 680000 \
    --hysteresis-k 5.0 \
    --t-end 1814400 \
    --orbit-semi-major-axis 6728000 \
    --orbit-inclination 1.13 \
    --angular-velocity 0.1,0.1,0.1 \
    --checkpoint-interval 450
