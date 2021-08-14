function hi()
  print("hi")
end

# The MIT License (MIT)
#
# Copyright (c) January, 2014 michael otte
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# a Jlist node (note that key is unused for Jlist opperations, but it is often
# helpful to have a key value associated with data)
mutable struct JListNode{T}
  child::JListNode{T}
  parent::JListNode{T}
  data::T
  key::Float64

  # constructors
  JListNode{T}() where {T} = new{T}()
end

# a simple Jlist
mutable struct JList{T}
  front::JListNode{T}
  back::JListNode{T}
  bound::JListNode{T} # bounds either side of the list
  length::Int

  # constructors
  function JList{T}() where {T}
    endNode::JListNode{T} = JListNode{T}()
    endNode.child = endNode
    endNode.parent = endNode
    new{T}(endNode, endNode, endNode, 0)
  end
end

function JlistPush(Jlist::TS, data::T) where {TS, T}



  newNode::JListNode{T} = JListNode{T}()

  newNode.parent = Jlist.front.parent
  newNode.child = Jlist.front

  if Jlist.length == 0
    Jlist.back = newNode
  else
    Jlist.front.parent = newNode
  end


  newNode.data = data
  Jlist.front = newNode
  Jlist.length +=1

#  data.ownJListNode = newNode # added for easy removal in center of Jlist
end

function JlistPush(Jlist::TS, data::T, key::Float64) where {TS, T}

  newNode = JListNode{T}()
  newNode.parent = Jlist.front.parent
  newNode.child = Jlist.front

  if Jlist.length == 0
    Jlist.back = newNode
  else
    Jlist.front.parent = newNode
  end

  newNode.data = data
  newNode.key = key
  Jlist.front = newNode
  Jlist.length +=1

#  data.ownJListNode = newNode # added for easy removal in center of Jlist
end

function JlistTop(Jlist::TS) where {TS}
  if Jlist.length == 0
    # Jlist is empty
    return false
  end
  return Jlist.front.data
end

function JlistTopKey(Jlist::TS) where {TS}
  if Jlist.length == 0
    # Jlist is empty
    return false
  end
  return (Jlist.front.data, Jlist.front.key)
end

function JlistPop(Jlist::TS) where {TS}
  if Jlist.length == 0
    # Jlist is empty
    return false
  end

  oldTop = Jlist.front
  if Jlist.length > 1
    Jlist.front.child.parent = Jlist.front.parent
    Jlist.front = Jlist.front.child
  elseif Jlist.length == 1
    Jlist.back = Jlist.bound
    Jlist.front = Jlist.bound
  end

  Jlist.length -=1

  oldTop.child = oldTop # added in case Jlist nodes hang around after this
  oldTop.parent = oldTop

  return oldTop.data
end

function JlistPopKey(Jlist::TS) where {TS}
  if Jlist.length == 0
    # Jlist is empty
    return false
  end

  oldTop = Jlist.front
  if Jlist.length > 1
    Jlist.front.child.parent = Jlist.front.parent
    Jlist.front = Jlist.front.child
  elseif Jlist.length == 1
    Jlist.back = Jlist.bound
    Jlist.front = Jlist.bound
  end

  Jlist.length -=1

  oldTop.child = oldTop # added in case Jlist nodes hang around after this
  oldTop.parent = oldTop
  return (oldTop.data, oldTop.key)
end


# removes thisNode from the list
function JlistRemove(Jlist::TS, thisNode::JListNode) where {TS}
  if Jlist.length == 0
    return true
  end


  if Jlist.front == thisNode
    Jlist.front = thisNode.child
  end
  if Jlist.back == thisNode
    Jlist.back = thisNode.parent
  end

  nextNode = thisNode.child
  previousNode = thisNode.parent

  if Jlist.length > 1 && previousNode != previousNode.child
    previousNode.child = nextNode
  end

  if Jlist.length > 1 && nextNode != nextNode.parent
   nextNode.parent = previousNode
  end

  Jlist.length -= 1

  if Jlist.length == 0
    Jlist.back = Jlist.bound # dummy node
    Jlist.front = Jlist.bound # dummy node
  end

  thisNode.parent = thisNode
  thisNode.child = thisNode

  return true
end



# assumes that the data has a subtype called ownJListNode that
# points to the JListNode that contains the data
# ALSO ASSUMES that there are not other external references
# to the JListNode that contains the data (it may break these if so)
#function JlistRemoveData{TS,T}(Jlist::TS, dataToRemove::T)
#  JlistRemove(Jlist, dataToRemove.ownJListNode)
#end




function JlistPrint(Jlist::TS) where {TS}
  ptr = Jlist.front
  while ptr != ptr.child
    println(ptr.data)
    ptr = ptr.child
  end
end

function JlistPrintKeys(Jlist::TS) where {TS}
  ptr = Jlist.front
  while ptr != ptr.child
    println(ptr.data)
    println("--------------: $(ptr.key)")
    ptr = ptr.child
  end
end

function JlistEmpty(Jlist::TS) where {TS}
  while JlistPop(Jlist) != false
  end
end

function testCase()
  L = JList{Array{Int64,2}}()

  a = [1 1 1]
  b = [2 2 2]
  c = [3 3 3]
  d = [4 4 4]

  JlistPush(L,a)
  JlistPush(L,b)
  JlistPush(L,c)
  JlistPrint(L)
  JlistEmpty(L)

  println("-")

  JlistPrint(L)

  println("-")

  JlistPush(L,a)
  JlistPush(L,b)
  JlistPush(L,c)
  JlistPush(L,b)
  JlistPush(L,a)
  cc = L.front
  JlistPrint(L)

  println("-")

  JlistRemove(L, cc)
  JlistPush(L,d)
  JlistPrint(L)
end
