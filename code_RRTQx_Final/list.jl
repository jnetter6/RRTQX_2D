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


# a list node (note that key is unused for list opperations, but it is often
# helpful to have a key value associated with data)
mutable struct ListNode{T}
  child::ListNode{T}
  data::T
  key::Float64

  # constructors
  ListNode{T}() where {T} = new{T}()
end

# a simple list
mutable struct List{T}
  front::ListNode{T}
  length::Int

  # constructors
  function List{T}() where {T}
    endNode::ListNode = ListNode{T}()
    endNode.child = endNode
    new{T}(endNode, 0)
  end
end

function listPush(list::TS, data::T) where {TS, T}
  newNode = ListNode{T}()
  newNode.child = list.front
  newNode.data = data
  list.front = newNode
  list.length +=1
end

function listPush(list::TS, data::T, key::Float64) where {TS, T}
  newNode = ListNode{T}()
  newNode.child = list.front
  newNode.data = data
  newNode.key = key
  list.front = newNode
  list.length +=1
end

function listTop(list::TS) where {TS}
  if list.front == list.front.child
    # list is empty
    return false
  end
  return list.front.data
end

function listTopKey(list::TS) where {TS}
  if list.front == list.front.child
    # list is empty
    return false
  end
  return (list.front.data, list.front.key)
end

function listPop(list::TS) where {TS}
  if list.front == list.front.child
    # list is empty
    return false
  end
  oldTop = list.front
  list.front = list.front.child
  list.length -=1
  return oldTop.data
end

function listPopKey(list::TS) where {TS}
  if list.front == list.front.child
    # list is empty
    return false
  end
  oldTop = list.front
  list.front = list.front.child
  list.length -=1
  return (oldTop.data, oldTop.key)
end


function listPrint(list::TS) where {TS}
  ptr = list.front
  while ptr != ptr.child
    println(ptr.data)
    ptr = ptr.child
  end
end

function listEmpty(list::TS) where {TS}
  while listPop(list) != false
  end
end

# helper used below
function listCopyGuts(list::TL, nodeExample::TN) where {TL, TN}
  newList = TL()

  ptr = list.front
  newList.front = TN()
  new_ptr = newList.front

  while ptr != ptr.child
    new_ptr.child = TN()
    new_ptr.data = ptr.data
    new_ptr.key = ptr.key

    new_ptr = new_ptr.child
    ptr = ptr.child
  end
  new_ptr.child = new_ptr
  #new_ptr.data = ptr.data
  new_ptr.key = ptr.key

  newList.length = list.length

  return newList
end

listCopy(list::TL) where {TL} = listCopyGuts(list, list.front)


function testCase()
  L = List{Array{Int64,2}}()

  a = [1 1 1]
  b = [2 2 2]
  c = [3 3 3]

  listPush(L,a)
  listPush(L,b)
  listPush(L,c)

  listPrint(L)
  listEmpty(L)
  listPrint(L)
  listPush(L,a)
  listPush(L,b)
  listPrint(L)

  println("-- copy:")
  L2 = listCopy(L)
  #listPush(L2,a)
  listPrint(L2)


  println("-- original:")
  #listPush(L,b)
  listPrint(L)

  println("$(L.length) =?= $(L2.length)")
end

function test()
  aa = 1
  bb = 2
  cc = 3
  L = List{Int64}()
  LL = List{Int64}()
  listPush(L, aa)
  listPush(L, bb)
  listPush(L, cc)
  listItem = L.front
  while (listItem != listItem.child)
    ob = deepcopy(listItem)
    listPush(L, ob.data)
    listPush(LL, ob.data)
    listItem = listItem.child
  end
  listPrint(L)
  println("---")
  listPrint(LL)
end
