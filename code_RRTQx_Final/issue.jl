

type Container
  data::Float64
end


function funC(n::Int)
  tic()
  C = Array(Container,n)
  for i = 1:n
    C[i] = Container(0.0)
  end
  println("  C's length: $(length(C))")
  println(" inside time: $(toq())")
end

function testC(n::Int = 100000000, manualGC::Bool = false, manualOffGC = false)
  if manualOffGC
    gc_disable()
  end
  tic()
  funC(n)
  println("outside time: $(toq())")
  if manualOffGC
    gc_enable()
  end
  if manualGC
    tic()
    gc()
    println(" manual gc(): $(toq()) (only time for this call)")
  end
  println("----------------------")
end




#
#
#julia> gc()
#
#julia> @time testC(10000000, false)
#  C's length: 10000000
# inside time: 0.263652254
#outside time: 0.263681618
#----------------------
#elapsed time: 0.26371345 seconds (240002000 bytes allocated, 48.61% gc time)
#
#julia> @time testC(10000000, true)
#  C's length: 10000000
# inside time: 0.258428296
#outside time: 0.258460679
# manual gc(): 0.047133994 (only time for this call)
#----------------------
#elapsed time: 0.30575224 seconds (240002728 bytes allocated, 56.58% gc time)



