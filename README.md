# Torso Tools

A set of helper functions to generate source and sensor layouts for magnetospinography simulations.

## Example usage

An example script can be found in `example/example.m` but perhaps its best to consult [this script](https://github.com/georgeoneill/study-spinevol/blob/main/sv_generate_geometries.m) to generate a source model and [this script](https://github.com/georgeoneill/study-spinevol/blob/main/sv_make_lead_fields_central.m) to generate the lead fields using a variety of forward models.

## Dependencies

torso_tools is dependent on [SPM](https://github.com/spm/spm) (not included) and an optional submodule [hbf_lc_p](https://github.com/MattiStenroos/hbf_lc_p) for boundary element analysis.

The repository also makes use of some modified assets from [ECGSim](https://www.ecgsim.org/downloads/other13/downloadgeo.php) which are shared under a GNU GPL v3 License.