using DelimitedFiles, CSV, DataFrames, PyCall, Distances, Plots, LibGEOS

function plot_pareto_front(objectives)
    p = Plots.scatter(objectives[:,1],objectives[:,2], xlabel="Coverage", ylabel="Preferences", title="Pareto front", legend=false)
    n = size(objectives,1)
    for i in 1:n
        annotate!(objectives[i,1], objectives[i,2], text(string(i), 8, :black, :bottom))
    end
    display(p)
end

function create_map(df, centers, π, radius)

    folium = pyimport("folium")
    # Legend = pyimport("from folium.plugins import Legend")
    
    map = df.explore(tiles="CartoDB positron",column=df["POPULATION"], cmap="magma_r")
    
    # Insert centroids
    n = df.shape[1]
    for i in 0:n-1
        lat = df.iloc[i].LATITUDE
        lon = df.iloc[i].LONGITUDE
        folium.Circle([lat, lon], 
                      weight = 1,
                      color="#000000",
                      radius=1).add_to(map)
    end
    
    # Insert centers
    for center in centers
        lat = df.loc[center].LATITUDE
        lon = df.loc[center].LONGITUDE
        if center in π
            folium.Circle([lat, lon], 
            radius=radius-245).add_to(map)
            #radius=radius).add_to(map)
        else
            folium.Circle([lat, lon], 
            color="#000000",
            radius=radius-245).add_to(map)
            #radius=radius).add_to(map) 
        end
    end

    # legend_colors = ["#000000","00FF00"]
    # legend_labels = ["Categoria 1","Categoria 2"]
    # legend = Legend(position="topright", colors=legend_colors, labels=legend_labels)
    # legend.add_to(map)
    
    return map
end

function load_geo_dataframe(geo_dataset)
    pd = pyimport("pandas")
    gpd = pyimport("geopandas")
    folium = pyimport("folium")

    df = gpd.read_file(geo_dataset)
    df["POPULATION"] = pd.to_numeric(df["POPULATION"])
    df["CENSUS_BLOCK_GROUP"] = gpd.GeoSeries.from_wkt(df["CENSUS_BLOCK_GROUP"])               
    df.set_geometry(col=df["CENSUS_BLOCK_GROUP"],inplace=true,crs="EPSG:4326")
    
    df["LATITUDE"] = pd.to_numeric(df["LATITUDE"])
    df["LONGITUDE"] = pd.to_numeric(df["LONGITUDE"])
    return df
end

function calculate_jaccard(C)
    return 1 .- pairwise(jaccard,C,dims=1)
end

function load_data_new(state, county, p, radius)
    
    #total_time = time()
    # Read distance matrix
    D = readdlm("data/$state" * "_" * "$county" * "_distance.csv")
    # Load geo dataframe
    df = load_geo_dataframe("data/$state" * "_" * "$county" * "_geo.csv")

    # Create cover matrix
    C = D .<= radius
    
    # Sort by population size and get ids permutation
    population = df["POPULATION"]
    id = sortperm(population.values,rev=true)
    #best_center_ids = id[1:p]
    
    #total_time = time() - total_time
    # println("Time to load dataset: $total_time")
        
    solucoes = execute_model_new(population.values, C, p)

    save_file = string(state,"_",county,"_",p,"_",radius,".txt")
    open(save_file, "w") do file
        id = 1
        for solucao in solucoes
            write(file, "###########\n")
            write(file, "Solution: $id\n")
            write(file, "Time: $(solucao.tempo)\n")
            write(file, "Objectives: $(solucao.objetivos)\n")
            write(file, "Centers: $(solucao.centers)\n")
            write(file, "Points: $(solucao.selected)\n")
            id += 1
        end
    end

    return solucoes
end



function load_data(state, county, p, radius)
    total_time = time()
    # Read distance matrix
    D = readdlm("data/$state" * "_" * "$county" * "_distance.csv")
    # Load geo dataframe
    df = load_geo_dataframe("data/$state" * "_" * "$county" * "_geo.csv")

    # Read distance matrix
    # D = readdlm("data_temp/$state" * "_" * "$county" * "_distance1.csv")
    # Load geo dataframe
    # df = load_geo_dataframe("data_temp/$state" * "_" * "$county" * "_geo1.csv")
    

    # LATITUDE = 34.11 34.01 -118.29 -118.19

    #map = df.explore(tiles="CartoDB positron",column=df["POPULATION"], cmap="magma_r")
    #map.save("saida.html")

    # Create cover matrix
    C = D .<= radius
    # Create Jaccard matrix
    # J = calculate_jaccard(C)
    # Sort by population size and get ids permutation
    population = df["POPULATION"]
    π = sortperm(population.values,rev=true)
    # Create precedence vector 'v' (each individual victory over others)
    n = size(π,1)
    v = zeros(Int64,n)
    for i in 1:n
        v[π[i]] = n-i 
    end
    # Create Borda Matrix
    BordaSum = zeros((n,n))
    for i in 1:n
        for j in 1:n
            BordaSum[i,j] = v[i] + v[j] 
        end
    end
    total_time = time() - total_time
    println("Time to load dataset: $total_time")

    #p = 125
    #p = 5
    #p = 50
    #p=100

    total_time = time()
    execute_model(BordaSum,C,p,π)
    total_time = time() - total_time
    println("Time to execute model: $total_time")

    return df, D, BordaSum, C, π
end

function calculate_intersection(points)
    a = area(points[1])
    n = size(points,1)
    I = zeros(n,n)
    for (i, g1) in enumerate(points)
        for (j, g2) in enumerate(points)
            # avoid counting double, and only use non empty intersections
            if i < j && LibGEOS.intersects(g1, g2)
                isect = LibGEOS.intersection(g1, g2)
                if typeof(isect) == Polygon
                    I[i,j] = I[j,i] = 100*(area(isect)/a)
                end
            end
        end
    end
    I = 100 .- I
    return I
end

function plot_solution_LibGEOS(circles, centers, best_ids)
    n = size(circles,1)
    plt = Plots.plot(circles, aspect_ratio=:equal,color=:blue, fillalpha=0.0)
    for i in 1:n
        if i in centers
            if i in best_ids
                Plots.plot!(circles[i], aspect_ratio=:equal,color=:red)
            else
                Plots.plot!(circles[i], aspect_ratio=:equal,color=:blue)
            end
        end
        center = LibGEOS.GeoInterface.coordinates(centroid(circles[i]))
        annotate!(center[1], center[2], text(string(i), 8, :black, :center))
    end
    display(plt)
end

function load_data_LibGEOS(state, county, p, radius)
    # Read distance matrix
    D = readdlm("data/$state" * "_" * "$county" * "_distance.csv")
    # Load geo dataframe
    df = load_geo_dataframe("data/$state" * "_" * "$county" * "_geo.csv")
    # Create cover matrix
    C = D .<= radius
    # Population
    population = df["POPULATION"]
    population = population.values

    # Define points in LibGEOS
    y = df["LONGITUDE"]
    x = df["LATITUDE"]
    n = size(C,1)
    points = [buffer(Point(x[i],y[i]), radius) for i in 1:n]

    total_time = time()
    # Find intersections between points radius
    I = calculate_intersection(points)
    total_time = time() - total_time
    #println(total_time)

    # Execute model
    execute_model_LibGEOS(points, population, C, I, p)

end