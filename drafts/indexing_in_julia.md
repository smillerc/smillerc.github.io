@def title = "Array indexing tips in Julia"
@def hascode = true
@def date = "5/18/2022"
@def tags = ["julia"]

# Array indexing tips in Julia

https://discourse.julialang.org/t/discussion-on-why-i-no-longer-recommend-julia-by-yuri-vishnevsky/81151/88

These are a collection of useful methods of iterating through arrays in Julia. Much of this is already in the manual or other documentation online, but this is sort of a cheatsheet for my own purposes. Hopefully this helps you.

One of the confusing aspects for me when I was starting in Julia how it handles arbitrary array interation. Coming from Fortran-land I was pleasantly surprised to see 


```
do j = 1, N
    do i = 1, M
        A(i,j) = A(i,j) + 1
    end do
end do
```


Using `eachindex`


CartesianIndices
```julia
for I in CartesianIndices(A)
    i,j = Tuple(I)
    @show i, j
end
```

OffsetArrays


```julia:diag1
using OffsetArrays

B = OffsetArray(reshape(1:9, 3, 3), -3, -2)
```

\show{diag1}

Loop along the diagonal of the `OffsetArray`
```julia:diag2
for (i, j) in zip(axes(B)...)
    @show i, j
end
```

\show{diag2}


