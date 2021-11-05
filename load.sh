#!/bin/bash
root -l <<EOF
.L ParticleType.cpp+
.L ResonanceType.cpp+
.L Particle.cpp+
.L GenerateParticles.cpp+
GenerateParticles();
.! cp histograms.root histograms_copy.root
.L AnalyzeData.cpp+
AnalyzeData();
EOF
