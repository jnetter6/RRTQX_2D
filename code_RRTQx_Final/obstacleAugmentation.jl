### this contains functions for obstacle augmentation
using LazySets
using DelimitedFiles
using Polyhedra
using DataFrames
# import LazySets.Approximations:overapproximate, approximate

function staticObstacleAugmentation(filename, augDist)
    # load obstacles
    fptr = open(filename, "r")

    # get the number of polygons
    P = parse(Int64, readline(fptr))
    # obstacle = Array{Array{Float64,2}, 1}(undef, P)
    array_V_obstacle = Array{Array{Array{Float64, 1}}, 1}(undef, P)

    for p = 1:P
        # get the number of point in this polygon
        N = parse(Int64, readline(fptr))

        # polygon = Array{Float64}(undef, N, 2)
        polygon = Array{Array{Float64, 1}}(undef, N)
        for n = 1:N
          polygon[n] = vec(str2array(readline(fptr)))
        end
        array_V_obstacle[p] = polygon
        # for i = 1:obsMult
        #   addObsToCSpace(S, Obstacle(3, polygon))
        # end
        readline(fptr)
    end
    close(fptr)

    # augment obstacles
    # array_V_obstacle = VPolygon.(array_obstacle)
    array_Vrep_obstacle = vrep.(array_V_obstacle) # Polyhedra.jl

    # V_obstacle1 = Array{VPolygon{Float64,Array{Float64,1}}, 1}(undef, P)
    # for p = 1:P V_obstacle1[p] = VPolygon(obstacle[p]) end
    # array_Ball_obstacle = LazySets.Approximations.ballinf_approximation.(array_V_obstacle)
    # array_Ball_obstacle1 = Array{BallInf{Float64}, 1}(undef, P)
    # for p = 1:P array_Ball_obstacle1[p] = LazySets.Approximations.ballinf_approximation(array_V_Obstacle[p]) end

    Ball_aug = Ball2(zeros(2), augDist)
    VPolygon_aug = convert(VPolygon, LazySets.Approximations.overapproximate(Ball_aug, .001))
    Vrep_aug = vrep(VPolygon_aug.vertices) # Polyhedra.jl
    array_allV_augmentedObstacle = array_Vrep_obstacle .+ Vrep_aug
    array_V_augmentedObstacle = convex_hull.(getfield.(array_allV_augmentedObstacle, :points))

    # write into new obstacle file





    # plot.(getindex.(array_obstacle, 1), getindex.(array_obstacle, 2))
    # p2 = plot!.(getindex.(array_Points_augmentedObstacle, 1), getindex.(array_Points_augmentedObstacle, 2))
    # array_V_aug = Array{VPolygon{Float64,Array{Float64,1}}, 1}(undef, P)
    # for p = 1:P array_V_aug[p] = V_aug end
    # start using Polyhedra.jl
    # array_V_augmentedObstacle = array_V_obstacle .⊕ array_V_aug
    # use LazySets.jl again to
    # array_V_augmentedObstacle = array_V_obstacle .⊕ array_V_aug
    # convert.(HPolygon, array_V_augmentedObstacle)
    # array_Ball_aug = Array{Ball2{Float64}, 1}(undef, P)
    # for p = 1:P array_Ball_aug[p] = Ball_aug end
    # # Ball_augmentedObstacle .= Ball_obstacle .+ Ball_aug
    # array_Ball_augmentedObstacle = array_Ball_obstacle .⊕ array_Ball_aug
    # array_H_augmentedObstacle = overapproximate.(array_Ball_augmentedObstacle, .01)
    # array_V_augmentedObstacle = convert.(VPolygon, array_H_augmentedObstacle)
    # # array_Vertices = Array{Array{Array{Float64, 1}, 1}, 1}(undef, P)
    # # for p = 1:P array_Vertices[p] = array_V_augmentedObstacle[p].vertices end
    # array_Vertices = getfield.(array_V_augmentedObstacle,:vertices)
    # array_v1 = Array{Array{Float64, 1}, 1}(undef, P)
    # array_v2 = Array{Array{Float64, 1}, 1}(undef, P)
    # for p = 1:P array_v1[p] = getindex.(array_Vertices[p], 1) end
    # for p = 1:P array_v2[p] = getindex.(array_Vertices[p], 2) end
    #
    # println(length(array_v1))
    # plot(array_H_augmentedObstacle[2], alpha=.2,aspectratio=1)


end

# save augmented polygon into ob.polygon, not change ob.polygonWithoutAug
function obstacleAugmentation(S::CSpace, augDist::Float64)
    Ball_aug = Ball2(zeros(2), augDist)
    VPolygon_aug = convert(VPolygon, LazySets.Approximations.overapproximate(Ball_aug, .01))
    Vrep_aug = vrep(VPolygon_aug.vertices) # Polyhedra.jl
    listItem = S.obstacles.front
    while (listItem != listItem.child)
        # if !listItem.data.obstacleUnused
        # obs = deepcopy(listItem)
        Vrep_obstacle = vrep(listItem.data.polygonWithoutAug)
        V_augmentedObstacle = Vrep_obstacle + Vrep_aug
        listItem.data.polygon = vcat(convex_hull(eachrow(V_augmentedObstacle.V)|>collect)'...)
        listItem.data.radius = sqrt(maximum(sum((broadcast(-, listItem.data.polygon, listItem.data.position)).^2, dims=2)))
        # listPush(S.obstacles, obs.data)
        # listPush(S.augObs, obs.data)
        # end
        listItem = listItem.child
    end
end

function showAugObs(oblist::List{Obstacle})
    pl = plot([-50,50,50,-50,-50], [-50,-50,50,50,-50], aspectratio = 1)
    listItem = oblist.front
    while (listItem != listItem.child)
        pl = [pl; plot!(listItem.data.polygon[:,1],listItem.data.polygon[:,2])]
        listItem = listItem.child
    end
    plot(pl...)
end

function obs()
    a = Ball2(zeros(2), 1.)
    # shape = [[-11.2788, -35.2193], [20.4032, -22.3538], [8.49078, -8.84503], [-9.75807, -18.1725]]
    shape = [[0,0],[0,1],[1,0]]
    p = LazySets.Approximations.ballinf_approximation(VPolygon([[-11.2788, -35.2193], [20.4032, -22.3538], [8.49078, -8.84503], [-9.75807, -18.1725]]))
    pshape = plot(getindex.(shape,1),getindex.(shape,2),aspectratio=1)
    ap = p⊕a

    inv = LazySets.Approximations.overapproximate(ap,.000001)
    # pl = plot(inv, alpha=0.2, aspectratio=1)
    vinv=convert(VPolygon,inv)
    v=vinv.vertices
    v1=getindex.(v,1)
    v2=getindex.(v,2)
    println(length(v1))
    paug = plot!(v1,v2,aspectratio=1)
    pbound = plot!([-50,50,50,-50],[-50,-50,50,50],aspectratio=1)
    plot(pshape,pbound,paug)
    # str = Array{String,1}(undef, 20)
    # for i=1:20 str[i]="$i" end
    # scatter!(v1,v2,series_annotations=str)
    # plot(p, aspectratio=1)
    # plot!(ap, alpha=0.6, aspectratio=1)
    # p0 = plot(b, 1e-6, aspectratio=1)
    # p1 = plot!(p0, overapproximate(b, 1.), alpha=0.4, aspectratio=1)
    #
    # p0 = plot(b, 1e-6, aspectratio=1)
    # p2 = plot!(p0, overapproximate(b, 0.1), alpha=0.4, aspectratio=1)
    #
    # p0 = plot(b, 1e-6, aspectratio=1)
    # p3 = plot!(p0, overapproximate(b, 0.01), alpha=0.4, aspectratio=1)
    #
    # plot(p1, p2, p3, layout=(1, 3))
    # plot!(Singleton(an_element(ap)))

    # println(getindex.(vinv.vertices,1)[1], getindex.(vinv.vertices,2)[1])
    # print(getindex.(vinv.vertices,1))
    # print(getindex.(vinv.vertices,2))
    # plot!.(getindex.(vinv.vertices,1), getindex.(vinv.vertices,2),"o")
end

function tb()
    b1 = LazySets.Approximations.overapproximate(Ball2(zeros(2), .57), .001)
    vb1 = convert(VPolygon, b1)
    # cb1 = convex_hull(vb1.vertices[:])
    # b2 = Ball2(zeros(2), .4)
    # b2 = VPolygon([[0.,0.],[0.,1.],[1.,0.]])
    b2 = VPolygon([[-51.324885, 34.897661], [-38.145161, 34.254386], [-39.665899, 28.786550], [-51.071429, 29.751462]])

    # b2 = convex_hull([[0.0,0.],[0.,1.],[1.,0.]])
    pb2 = plot(b2, aspectratio=1)
    b = vb1⊕b2
    pb=plot!(b, aspectratio=1)
    plot(pb2,pb)
    # bb = overapproximate(b, .1)
    # plot!(bb, aspectratio=1)
end

function tt()
    aa = [1 2; 3 4; 5 6]
    bb = Array{Array{Float64, 1}}(undef, size(aa, 1))
    for i = 1:size(aa, 1) bb[i] = aa[i, :] end
    println(bb)
end
