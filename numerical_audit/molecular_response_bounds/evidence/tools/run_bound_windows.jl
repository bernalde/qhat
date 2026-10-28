# Root-owned two-case dispatcher; no fresh enumeration and no cap reset.
using TOML, SHA
length(ARGS)==2 || error("usage: run_bound_windows.jl BOUNDS_CODE NEW_OUTPUT_DIRECTORY")
code=abspath(ARGS[1]);output=abspath(ARGS[2])
ispath(output) && error("fresh output directory required")
include(joinpath(code,"window_development.jl"))
mkpath(output)
archive=joinpath(code,"..","julia_cubic_feasibility")
for (short,case,start) in (("CH4","CH4_r1p09",157),("CO","CO_r1p13",168))
    input=joinpath(archive,"evidence","inputs","cubic",case*".toml")
    result=joinpath(output,short*".BOUND.a1.toml")
    try
        record=development(input,start,result)
        oraclepath=joinpath(archive,"evidence","cubic_dp_validation","results",short*"_MAIN.a1.toml")
        oracle=TOML.parsefile(oraclepath)
        oracle["status"]=="complete" && oracle["input_sha256"]==filehash(input) || error("archived input/oracle identity differs")
        oracle["flat_sha256"]==record["flat_sha256"] || error("frozen table bits differ")
        oracle["window_positions"]==collect(start:start+7) || error("oracle window differs")
        r=Int(record["r_inc"])
        r==oracle["steps"] || error("oracle affordability differs")
        optimum=rationalparse(oracle["expected_Q"])
        0<=record["L_Q"]<=optimum<=record["Q_inc"] || error("bound/oracle ordering differs")
        optimum/r^4==rationalparse(oracle["expected_J"]) || error("oracle J normalization differs")
        write_record(joinpath(output,short*".ORACLE_CHECK.a1.toml"),Dict(
            "status"=>"verified_restricted_bound_against_archived_exhaustive_optimum",
            "input_sha256"=>filehash(input),"bound_sha256"=>filehash(result),
            "oracle_sha256"=>filehash(oraclepath),"flat_sha256"=>record["flat_sha256"],
            "L_Q"=>record["L_Q"],"optimum_Q"=>optimum,"incumbent_Q"=>record["Q_inc"],
            "steps"=>r,"full_molecule_ratio_claimed"=>false,"fresh_enumeration"=>false))
        println("WINDOW_BOUND_CHECKED ",case);flush(stdout)
    catch exception
        write_record(joinpath(output,short*".REFUSAL.a1.toml"),Dict(
            "status"=>"incomplete_window_bound_or_oracle_check","reason"=>sprint(showerror,exception),
            "input_sha256"=>filehash(input),"positive_full_bound"=>false,
            "bound_artifact_retained"=>isfile(result)))
        println("WINDOW_BOUND_INCOMPLETE ",case," ",typeof(exception));flush(stdout)
    end
end
