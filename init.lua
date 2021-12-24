--[[

	Orkan
	@myteawastaken
	---
	Orkan._g2p2g_2(domain, p0, p1)
	Orkan._g2p2g_3(domain, p0, p1)
	-
	Orkan.advanceParticles(domain, p0, p1)
	Orkan.advanceGrid(domain)
	Orkan.createDomain(dt, size, constantEulerianExternalForces)
	Orkan.addParticle(domain, position, mass, material)
]]

local Math = require(script.Math)

local clamp = math.clamp
local v2, v3 = Vector2.new, Vector3.new
local v2b, v3b = v2(), v3()

--
-- Generates 3^dim relative cell positions, placed in order in array {Vector(i, j[, k])}. Supports two and three dimensions.
--
local function generateNeighbouringCells(dim, min, max)
	local neighbouringCells = {}
	for i = min, max do
		for j = min, max do
			if dim == 3 then
				for k = min, max do
					table.insert(neighbouringCells, v3(i, j, k))
				end
				continue
			end
			table.insert(neighbouringCells, v2(i, j))
		end
	end
	return neighbouringCells
end

--
-- Improve algorithm performance with ipairs and already cached vectors
--
local neighbouringCells2 = generateNeighbouringCells(2, 0, 2)
local neighbouringCells3 = generateNeighbouringCells(3, 0, 2)

--
-- Get number of dimensions of vector
--
local function getDim(vec)
	return if typeof(vec) == 'Vector2' then 2 else 3
end

--
-- Get number of dimensions of vector
--
local function getVolume(vec)
	local volume = vec.X * vec.Y
	if typeof(vec) == 'Vector3' then
		volume *= vec.Z
	end
	return volume
end

--
-- Computes and applies weight gradient to array "weightGradient" from distance "distance" (2 dimensions)
--
local v2_1 = v2(1, 1)
local v2_1_5, v2_0_75, v2_0_5 = v2_1 * 1.5, v2_1 * 0.75, v2_1 * 0.5
local function computeWeightGradient2(weightGradient, distance)
	weightGradient[1] = 0.5 * Math.square(v2_1_5 - distance)
	weightGradient[2] = v2_0_75 - Math.square(distance - v2_1)
	weightGradient[3] = 0.5 * Math.square(distance - v2_0_5)
end

--
-- Computes and applies weight gradient to array "weightGradient" from distance "distance" (3 dimensions)
--
local v3_1 = v3(1, 1, 1)
local v3_1_5, v3_0_75, v3_0_5 = v3_1 * 1.5, v3_1 * 0.75, v3_1 * 0.5
local function computeWeightGradient3(weightGradient, distance)
	weightGradient[1] = 0.5 * Math.square(v3_1_5 - distance)
	weightGradient[2] = v3_0_75 - Math.square(distance - v3_1)
	weightGradient[3] = 0.5 * Math.square(distance - v3_0_5)
end

local Orkan = {}

local function vec2ToInt(vec, boundX)
	return vec.X + boundX * vec.Y
end

local function vec3ToInt(vec, boundX, boundY)
	return vec.X + boundX * vec.Y + boundX * boundY * vec.Z
end

--
-- Runs 2D G2P2G on particles [p0] to [p1] in domain "domain" with initialized weight gradient "weightGradient"
--
function Orkan._g2p2g_2(domain, p0, p1, weightGradient)
	local x, ms, mt, F = domain.x, domain.ms, domain.mt, domain.F
	local gmOut, gvOut, gvIn = domain.gmOut, domain.gvOut, domain.gvIn
	local size, dt = domain.size, domain.dt
	local nParticles = domain.nParticles
	local sx, sy = size.X, size.Y
	for p = p0, p1 do
		local xp = x[p]

		local baseCell, dx
		-- Thanks to the G2P2G pipeline, we can interpolate velocity for the grid, hence storing particle velocity is not needed
		local vp, Cp = v2b, Math.blankMatrix(2, 2)

		-- We don't want to waste a G2P for new particles (this is timestep k+1!)
		if p < nParticles then
			--
			-- G2P (timestep k+1)
			--
			baseCell = Math.floorVec2(xp - v2_0_5)
			dx = xp - baseCell

			computeWeightGradient2(weightGradient, dx)

			for _, offset in ipairs(neighbouringCells2) do
				-- TODO: make sure this is correct
				local dxi = offset - dx
				local i = vec2ToInt(baseCell + offset, sx) + 1
				local weight = weightGradient[offset.X + 1].X * weightGradient[offset.Y + 1].Y
				local vi = gvIn[i]
				vp += weight * vi
				Cp = Math.appliedMat22Add(Cp, Math.matScalMul(Math.outerProduct22(vi, dxi), 4 * weight))
			end
			xp += dt * vp
			xp = v2(
				clamp(xp.X, 1, sx - 1),
				clamp(xp.Y, 1, sy - 1)
			)
			x[p] = xp


		else
			nParticles += 1
		end
		--
		-- P2G (timestep k)
		--
		baseCell = Math.floorVec2(xp - v2_0_5)
		dx = xp - baseCell

		computeWeightGradient2(weightGradient, dx)

		local mp = ms[p]
		local affine = Cp

		for _, offset in ipairs(neighbouringCells2) do
			local dxi = offset - dx
			local i = vec2ToInt(baseCell + offset, sx) + 1
			local weight = weightGradient[offset.X + 1].X * weightGradient[offset.Y + 1].Y
			gvOut[i] += weight * (mp * vp + Math.matVecMul22(affine, dxi))
			gmOut[i] += weight * mp
		end
	end
	domain.nParticles = nParticles
end

--
-- Runs 2D G2P2G on particles [p0] to [p1] in domain "domain" with initialized weight gradient "w"
--
function Orkan._g2p2g_3(domain, p0, p1, w)
	for p = p0, p1 do
	end
end




--
-- Advances particles [p0] to [p1] in domain "domain" with initialized weight gradient "weightGradient"
--
function Orkan.advanceParticles(domain, p0, p1, weightGradient)
	if domain.dim == 2 then
		Orkan._g2p2g_2(domain, p0, p1, weightGradient)
	else
		Orkan._g2p2g_3(domain, p0, p1, weightGradient)
	end
end




--
-- Converts momentum to velocity, applies constant external eulerian forces and swaps grid input/output
--
function Orkan.advanceGrid(domain)
	local gmOut, gvOut, gvIn = domain.gmOut, domain.gvOut, domain.gvIn
	local extForces = domain.constantEulerianExternalForces
	local dt = domain.dt
	local blankVec = 0 * domain.size
	-- Velocity -> momentum and eulerian external forces
	for i, mass in ipairs(gmOut) do
		if mass > 0 then
			gvIn[i] = gvOut[i] / mass + dt * extForces
			gvOut[i] = blankVec
			gmOut[i] = 0
		end
	end
end




--
-- Create new orkan domain
--
function Orkan.createDomain(dt, size, constantEulerianExternalForces)
	local blank, n = size * 0, getVolume(size)
	return {
		---
		dim = getDim(size),		-- Domain's number of dimensinos
		size = size,			-- Domain's grid size
		---
		gmOut = table.create(n, 0),		-- Output mass tensor field (written to during each timestep, cleared between timesteps)
		gvOut = table.create(n, blank),	-- Output velocity tensor field (written to during each timestep, cleared between timesteps)
		gvIn = table.create(n, blank),	-- Input velocity tensor field (output velocity tensor field from last timestep)
		---
		x = {},		-- Particle positions
		ms = {},	-- Particle masses
		mt = {},	-- Particle materials
		F = {},		-- Particle F-matrices
		---
		nParticles = 0,		-- Number of particles (used and set internally and internally ONLY)
		---
		constantEulerianExternalForces = constantEulerianExternalForces,	-- The contant external forces applied to grid nodes.
		---
		dt = dt		-- Delta-time used for ex. advection (sim stability depends on this, lower dt=more stable)
		---
	}
end




--
-- Adds particle with attributes [position, mass, material] to orkan domain
--
function Orkan.addParticle(domain, position, mass, material)
	local n = #domain.x + 1
	domain.F[n] = Math.identityMatrix(domain.dim)
	domain.ms[n] = mass
	domain.mt[n] = material
	domain.x[n] = position
end




--
-- Bulk adds particles with attributes [positions[p], masses[p], materials[p]] to orkan domain
--
function Orkan.bulkAddParticles(domain, positions, masses, materials)
	local x, ms, mt, F = domain.x, domain.ms, domain.mt, domain.F
	local identityMat = Math.identityMatrix(domain.dim)
	local n = #x
	for p, xp in ipairs(positions) do
		n += 1
		F[n] = identityMat
		ms[n] = masses[p]
		mt[n] = materials[p]
		x[n] = xp
	end
end



return Orkan