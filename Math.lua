local floor = math.floor
local v2, v3 = Vector2.new, Vector3.new

local Math = {}

--
-- Removes decimals of vector elements of a 2-dimensional vector
--
function Math.floorVec2(vec)
	return v2(
		floor(vec.X),
		floor(vec.Y)
	)
end


--
-- Squares object/number "v". Works as long as object has metamethod "__mul" defined.
--
function Math.square(v)
	return v * v
end


--
-- Creates blank matrix of dimensions [i, j]
--
function Math.blankMatrix(i, j)
	return table.create(i * j, 0)
end
--
-- Creates identity matrix of dimensions [s, s]
--
function Math.identityMatrix(s)
	local matrix = table.create(s * s, 0)
	-- For each x, set matrix[x, y] (x==y) to 1
	for i = 0, s - 1 do
		matrix[i + s * i + 1] = 1
	end
	return matrix
end


--
-- Computes outer product of two 2-dimensional vectors v0 and v1
--
function Math.outerProduct22(v0, v1)
	return {
		v0.X*v1.X, v0.X*v1.Y,
		v0.Y*v1.X, v0.Y*v1.Y
	}
end


--
-- Computes addition of two 2x2 matrices m0 and m1
--
function Math.mat22Add(m0, m1)
	return {
		m0[1] + m1[1], m0[2] + m1[2],
		m0[3] + m1[3], m0[4] + m1[4]
	}
end
--
-- Computes addition of two 2x2 matrices m0 and m1. Writes result to matrix m0.
--
function Math.appliedMat22Add(m0, m1)
	for i in ipairs(m0) do
		m0[i] += m1[i]
	end
	return m0
end


--
-- Computes multiplication of 2x2 matrix "mat" and 2-dimensional vector "vec"
--
function Math.matVecMul22(mat, vec)
	return v2(
		mat[1]*vec.X + mat[2]*vec.X,
		mat[3]*vec.Y + mat[4]*vec.Y
	)
end


--
-- Computes multiplication of matrix "mat" and scalar "scal"
--
function Math.matScalMul(mat, scal)
	local res = {}
	for i, v in ipairs(mat) do
		res[i] = v * scal
	end
	return res
end
--
-- Computes multiplication of matrix "mat" and scalar "scal". Writes result to matrix "mat".
--
function Math.appliedMatScalMul(mat, scal)
	for i, v in ipairs(mat) do
		mat[i] *= scal
	end
end

return Math