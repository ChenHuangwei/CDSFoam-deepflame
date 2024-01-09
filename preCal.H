// --- Directed interpolation of primitive fields onto faces

volScalarField rPsi("rPsi", 1.0/psi);
surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
