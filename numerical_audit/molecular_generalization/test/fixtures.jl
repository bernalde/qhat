using SHA
fixture_digest(text)=bytes2hex(sha256(text))
function fixture_manifest(;mode="main",unit_count=2,n=12,electrons=8)
    reps=mode=="main" ? 2 : 1
    manifest=Dict{String,Any}("version"=>1,"order_index_base"=>1,"case"=>"fixture","n"=>n,"electrons"=>electrons,
        "partition"=>"conserving","input_sha256"=>fixture_digest("model"),"source_json_sha256"=>fixture_digest("json"),
        "protocol_sha256"=>fixture_digest("protocol"),"source_freeze_sha256"=>fixture_digest("freeze"),
        "mode"=>mode,"repetitions"=>reps,"expected_logical_rows"=>36reps,"training_time"=>1.0,
        "training_budget_multiplier"=>30,"times"=>[.65,1.35],"rotation_cap_multipliers"=>[8,24],"candidates"=>Any[])
    for rep in 1:reps, (arm,cpu) in (("parent",0),("canonical",0),("magnitude",0),("full",32),("full",64),
            ("window",32),("window",64),("annealing",32),("annealing",64)), t in (.65,1.35), cap in (8,24)
        oid="rep$(rep).$(arm).B$(cpu)"; digest=fixture_digest(oid)
        row=Dict{String,Any}("name"=>oid*".T$(t).R$(cap)","order_id"=>oid,"arm"=>arm,"repetition"=>rep,
            "cpu_budget"=>cpu,"partition"=>"conserving","formula"=>2,"time"=>t,"budget_multiplier"=>cap,
            "order"=>collect(1:unit_count),"constructor_sha256"=>digest,"claimed_score"=>1e-6,
            "score_kind"=>"diagonal_reference_cubic_response","eligible"=>true)
        arm in ("canonical","magnitude") || (row["parent_sha256"]=fixture_digest("rep$(rep).parent.B0"))
        push!(manifest["candidates"],row)
    end
    manifest
end
