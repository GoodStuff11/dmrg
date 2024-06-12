using YAML

"""
Reads commandline arguments (if any). After the file, the command line should
have the following form:

VARNAME value VARNAME2 value2 VARNAME3 value3

This will set all values from the YAML data_file (input_quick.yml) to the value specified
"""
function get_input_data(data_file::String="input_quick.yml"; default_filename::String)
    data = YAML.load_file(data_file)
    if length(ARGS) > 0
        for i in 1:2:length(ARGS)
            try
                data[ARGS[i]] = eval(Meta.parse(ARGS[i+1]))
            catch
                # in cases where the input is a string
                data[ARGS[i]] = ARGS[i+1]
            end
        end
    end

    mmax = data["m"]
    Nsites = data["Nsites"]
    Nbonds = data["Nbond"]
    Nsweep = data["Nsweep"]
    e_cutoff = data["ecutoff"]
    SVD_error = data["SVD"]
    gstart = data["gstart"]
    delta_g = data["dg"]
    Ng = data["Ng"]
    mbond = data["mbond"]
    pairs = data["pairs"]
    evod = data["states"]
    angle = data["angle"]*pi/180.0
    Estrength = data["strength"]
    Nstates = data["Nstates"]
    output_filename = get!(data, "filename", default_filename)
    parity_symmetry = get!(data, "ParitySymmetry", "even")
    inversion_symmetry = get!(data,"InversionSymmetry",  "none")

    return mmax, Nsites, Nbonds, Nsweep, e_cutoff, SVD_error, gstart, delta_g, Ng,
           mbond, pairs, evod, angle, Estrength, Nstates, output_filename, parity_symmetry,
           inversion_symmetry
end