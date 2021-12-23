local floor = math.floor
local v2, v3 = Vector2.new, Vector3.new

export type matrix = {number}

local Math = {}

function Math.SquareVector(vec: Vector2 | Vector3): Vector2 | Vector3
    return vec * vec
end

--- Converts vector2 "vec" to integer, given maximum X axis position "max_x" (= vector may only consist of positive integers)
function Math.BoundedVec2ToInt(vec: Vector2, max_x: number): number
    return vec.X + max_x * vec.Y
end

function Math.FloorVec2(vec: Vector2): Vector2
    return v2(
        floor(vec.X),
        floor(vec.Y)
    )
end

--- Creates a blank matrix of size [i, j]
function Math.BlankMatrix(i: number, j: number): matrix
    return table.create(i * j, 0)
end

--- Creates an identity matrix of size [s, s]
function Math.IdentityMatrix(s: number): matrix
local matrix = table.create(s * s, 0)
    for i = 0, s - 1 do
        matrix[i + s * i + 1] = 1
    end
    return matrix
end

--@source https://en.wikipedia.org/wiki/Outer_product
function Math.Mat22OuterProd(v0: Vector2, v1: Vector2): matrix
    return {
        v0.X*v1.X, v0.X*v1.Y,
        v0.Y*v1.X, v0.Y*v1.Y
    }
end

--- Computes addition of two 2x2 matrices m0 and m1
function Math.Mat22Add(m0: matrix, m1: matrix): matrix
    return {
        m0[1] + m1[1], m0[2] + m1[2],
        m0[3] + m1[3], m0[4] + m1[4]
    }
end

--- Computes addition of two 2x2 matrices m0 and m1, applies result to matrix m0
function Math.AppliedMat22Add(m0: matrix, m1: matrix): matrix
    for i in ipairs(m0) do
        m0[i] += m1[i]
    end
    return m0
end

--- Computes addition of two 3x3 matrices m0 and m1
function Math.Mat33Add(m0: matrix, m1: matrix): matrix
    return {
        m0[1] + m1[1], m0[2] + m1[2], m0[3] + m1[3],
        m0[4] + m1[4], m0[5] + m1[5], m0[6] + m1[6],
        m0[7] + m1[7], m0[8] + m1[8], m0[9] + m1[9]
    }
end

function Math.Mat22Vec2Mul(mat: matrix, vec: Vector2): Vector2
    return v2(
        mat[1]*vec.X + mat[2]*vec.Y,
        mat[3]*vec.X + mat[4]*vec.Y
    )
end

--- Computes product of a scalar and a matrix, where the matrix may be of any dimensions.
function Math.MatScalMul(scalar: number, mat: matrix): matrix
    local product = table.create(#mat) -- TODO: Should this change to just "{}"?
    for i, n in ipairs(mat) do
        product[i] = n * scalar
    end
    return product
end

--- Computes product of a scalar and a matrix, where the matrix may be of any dimensions. Avoids table creation by overriding the original matrix, but may therefore not be ideal in all cases.
function Math.AppliedMatScalMul(scalar: number, mat: matrix): matrix
    for i in ipairs(mat) do
        mat[i] *= scalar
    end
    return mat -- For convenience it still returns result, as it in certain cases may "look better".
end

return Math