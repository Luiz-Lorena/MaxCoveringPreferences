include("src/MaxCoverPreferences.jl")
using .MaxCoverPreferences

parameters = [
    "48" "113" 10;
    "06" "059" 10;
    "06" "073" 10;
    "17" "031" 10;
    "06" "037" 10;
    "NY" "NY" 10
]

for i in 1:size(parameters,1)
    println(parameters[i,:])
    MaxCoverPreferences.load_data_new(parameters[i,1],parameters[i,2],parameters[i,3],9000);
end