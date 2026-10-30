# OLD fixed-order mathematical bridge; no order optimization or physical errors.
using TOML,SHA
length(ARGS)==4 || error("OLD_KERNEL NEW_KERNEL ORIGINAL_REQUEST FRESH_OUTPUT")
const OLD_KERNEL=abspath(ARGS[1])
const NEW_KERNEL=abspath(ARGS[2])
const REQUEST=abspath(ARGS[3])
const OUTPUT=abspath(ARGS[4])
!ispath(OUTPUT)||error("fresh bridge output required")
mkdir(OUTPUT)
module OldBridge
include(Main.OLD_KERNEL)
end
module NewBridge
include(Main.NEW_KERNEL)
end
const O=OldBridge.FullOrderAnnealing
const N=NewBridge.FullOrderAnnealing
rational_text(x)=begin;a,b=split(x,'/');parse(BigInt,a)//parse(BigInt,b);end
request=TOML.parsefile(REQUEST)
score_close(a,b)=isfinite(a)&&isfinite(b)&&abs(a-b)<=max(1e-28,abs(a)*1e-8)
function units_for(A,records)
    UnitType=getfield(A,:Unit)
    result=Vector{UnitType}()
    for u in records
        exact=rational_text.(u["coefficients"])
        push!(result,UnitType(UInt64(u["flip"]),UInt64.(u["parities"]),Int8.(u["signs"]),Float64.(exact);exact=exact))
    end
    result
end
function save_report(path,report)
    !ispath(path)||error("preserve bridge artifact")
    open(path,"w") do io;TOML.print(io,report;sorted=true);end
end
report=Dict{String,Any}("schema"=>"molecular-old-scalar-cache-bridge-v1","complete"=>false,
    "old_source_sha256"=>bytes2hex(sha256(read(OLD_KERNEL))),"new_source_sha256"=>bytes2hex(sha256(read(NEW_KERNEL))),
    "request_sha256"=>bytes2hex(sha256(read(REQUEST))),"scope"=>"Two fixed old LiH orders; no search or physical errors",
    "rows"=>Dict{String,Any}[])
save_report(joinpath(OUTPUT,"BRIDGE.START.toml"),report)
old_context=O.FullOrderContext(units_for(O,request["units"]),request["n"],request["reference"],request["budget"])
new_context=N.FullOrderContext(units_for(N,request["units"]),request["n"],request["reference"],request["budget"])
orders=[("measured_parent",Int.(request["order"]).+1),("canonical",collect(1:length(old_context.units)))]
for (index,(label,order)) in enumerate(orders)
    before=O.full_column(old_context,order).column;after=N.full_column(new_context,order).column
    destinations=sort!(collect(union(keys(before),keys(after))))
    oldvalues=[get(before,state,0.)for state in destinations]
    newvalues=[get(after,state,0.)for state in destinations]
    bitwise=all(reinterpret(UInt64,a)==reinterpret(UInt64,b)for(a,b)in zip(oldvalues,newvalues))
    scale=isempty(oldvalues) ? 0. : maximum(abs,oldvalues)
    difference=isempty(oldvalues) ? 0. : maximum(abs.(oldvalues.-newvalues))
    column_verified=difference<=max(3e-14,scale*1e-10)
    oldscore=O.full_score!(old_context,order);newscore=N.full_score!(new_context,order)
    row=Dict{String,Any}("label"=>label,"order_zero_based"=>order.-1,"columns_bitwise_equal"=>bitwise,
        "columns_verified"=>column_verified,"column_max_abs_difference"=>difference,
        "old_score"=>Dict(string(k)=>v for(k,v)in pairs(oldscore)),"new_score"=>Dict(string(k)=>v for(k,v)in pairs(newscore)),
        "score_verified"=>oldscore.steps==newscore.steps&&score_close(oldscore.value,newscore.value)&&score_close(oldscore.cubic_response,newscore.cubic_response),
        "column_states"=>Int.(destinations),"old_column"=>oldvalues,"new_column"=>newvalues)
    push!(report["rows"],row)
    save_report(joinpath(OUTPUT,"BRIDGE.PARTIAL."*string(index)*".toml"),report)
    println("OLD_CACHE_BRIDGE_ROW ",label," bitwise ",bitwise," verified ",column_verified&&row["score_verified"]);flush(stdout)
end
report["old_kernel_stats"]=copy(old_context.stats);report["new_kernel_stats"]=copy(new_context.stats)
report["new_basis_ranks"]=length.(new_context.scalar_bases)
report["complete"]=all(row["columns_verified"]&&row["score_verified"]for row in report["rows"])
save_report(joinpath(OUTPUT,"BRIDGE.toml"),report)
exit(report["complete"] ? 0 : 2)
