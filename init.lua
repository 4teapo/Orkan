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

local neighbours2 = {
	v2(0, 0),
	v2(0, 1),
	v2(0, 2),
	v2(1, 0),
	v2(1, 1),
	v2(1, 2),
	v2(2, 0),
	v2(2, 1),
	v2(2, 2)
}

local neighbours3 = {
	v3(0, 0, 0),
	v3(0, 0, 1),
	v3(0, 0, 2),
	v3(0, 1, 0),
	v3(0, 1, 1),
	v3(0, 1, 2),
	v3(0, 2, 0),
	v3(0, 2, 1),
	v3(0, 2, 2),
	v3(1, 0, 0),
	v3(1, 0, 1),
	v3(1, 0, 2),
	v3(1, 1, 0),
	v3(1, 1, 1),
	v3(1, 1, 2),
	v3(1, 2, 0),
	v3(1, 2, 1),
	v3(1, 2, 2),
	v3(2, 0, 0),
	v3(2, 0, 1),
	v3(2, 0, 2),
	v3(2, 1, 0),
	v3(2, 1, 1),
	v3(2, 1, 2),
	v3(2, 2, 0),
	v3(2, 2, 1),
	v3(2, 2, 2)
}

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
	weightGradient[1] = 0.5 * Math.SquareVector(v2_1_5 - distance)
	weightGradient[2] = v2_0_75 - Math.SquareVector(distance - v2_1)
	weightGradient[3] = 0.5 * Math.SquareVector(distance - v2_0_5)
end

--
-- Computes and applies weight gradient to array "weightGradient" from distance "distance" (3 dimensions)
--
local v3_1 = v3(1, 1, 1)
local v3_1_5, v3_0_75, v3_0_5 = v3_1 * 1.5, v3_1 * 0.75, v3_1 * 0.5
local function computeWeightGradient3(weightGradient, distance)
	weightGradient[1] = 0.5 * Math.SquareVector(v3_1_5 - distance)
	weightGradient[2] = v3_0_75 - Math.SquareVector(distance - v3_1)
	weightGradient[3] = 0.5 * Math.SquareVector(distance - v3_0_5)
end

local Orkan = {}

--
-- Runs 2D G2P2G on particles [p0] to [p1] in domain "domain" with initialized weight gradient "w"
--
function Orkan._g2p2g_2(domain, p0, p1, weightGradient)
	local x, ms, mt, F = domain.x, domain.ms, domain.mt, domain.F
	local gm_out, gv_out, gv_in = domain.gm_out, domain.gv_out, domain.gv_in
	local size, dt = domain.size, domain.dt
	local nParticles = domain.nParticles
	local sx, sy = size.X, size.Y
	for p = p0, p1 do
		local xp = x[p]

		local baseCell, dx
		--
		-- Thanks to the G2P2G pipeline, we can interpolate velocity for the grid, hence storing particle velocity is not needed
		--
		local vp, Cp = v2b, Math.BlankMatrix(2, 2)

		--
		-- We don't want to waste a G2P for new particles (this is timestep k+1!)
		--
		if p < nParticles then
			--
			-- G2P (timestep k+1)
			--
			baseCell = Math.FloorVec2(xp)
			dx = xp - baseCell

			computeWeightGradient2(weightGradient, dx)

			for _, offset in ipairs(neighbours2) do
				-- TODO: make sure this is correct
				local dxi = offset - dx
				local i = Math.BoundedVec2ToInt(baseCell + offset, sx) + 1
				local w = weightGradient[offset.X + 1].X * weightGradient[offset.Y + 1].Y
				local vi = gv_in[i]
				vp += w * vi
				Cp = Math.AppliedMat22Add(Cp, Math.MatScalMul(4 * w, Math.Mat22OuterProd(vi, dxi)))
			end
		else
			nParticles += 1
		end
		--
		-- Advect particle
		--
		xp += dt * vp
		xp = v2(
			clamp(xp.X, 1, sx - 1),
			clamp(xp.Y, 1, sy - 1)
		)
		--
		-- P2G (timestep k)
		--
		baseCell = Math.FloorVec2(xp)
		local dx = xp - baseCell

		computeWeightGradient2(weightGradient, dx)

		local mp = ms[p]
		local affine = Cp

		for _, offset in ipairs(neighbours2) do
			local dxi = offset - dx
			local i = Math.BoundedVec2ToInt(baseCell + offset, sx) + 1
			local w = weightGradient[offset.X + 1].X * weightGradient[offset.Y + 1].Y
			gv_out[i] += w * (mp * vp + Math.Mat22Vec2Mul(affine, dxi))
			gm_out[i] += w * mp
		end
	end
end

function Orkan._g2p2g_3(domain, p0, p1, w)
	for p = p0, p1 do
	end
end

--
-- Advances particles [p0] to [p1] in domain "domain" with initialized weight gradient "w"
--
function Orkan.advanceParticles(domain, p0, p1, w)
	if domain.dim == 2 then
		Orkan._g2p2g_2(domain, p0, p1, w)
	else
		Orkan._g2p2g_3(domain, p0, p1, w)
	end
end

--
-- Converts momentum to velocity, applies constant external eulerian forces and swaps grid input/output
--
function Orkan.advanceGrid(domain)
	local gm_out, gv_out, gv_in = domain.gm_out, domain.gv_out, domain.gv_in
	local extForces = domain.constantEulerianExternalForces
	local dt = domain.dt
	-- Velocity -> momentum and eulerian external forces
	for i, mass in ipairs(gm_out) do
		if mass > 0 then
			local v = gv_in[i]
			gv_in[i] = gv_out[i] / mass + dt * extForces
			gv_out[i] = v
			gm_out[i] = 0
		end
	end
end


--
-- Create new orkan instance
--
function Orkan:createInstance(pip)
	return {
		domain = nil
	}
end


--
-- Create new orkan domain
--
function Orkan.createDomain(dt, size, constantEulerianExternalForces)
	local blank, n = size * 0, getVolume(size)
	return {
		---
		dim = getDim(size),
		size = size,
		---
		gmOut = table.create(n, 0),
		gvOut = table.create(n, blank),
		gvIn = table.create(n, blank),
		---
		x = {},
		ms = {},
		mt = {},
		F = {},
		---
		nParticles = 0,
		---
		constantEulerianExternalForces = constantEulerianExternalForces,
		dt = dt
		---
	}
end


--
-- Bind domain to orkan instance
--
function Orkan:bindDomainToInstance(instance, domain)
	instance.domain = domain
	instance.pip = self._pip(domain.dim)
end


--
-- Adds particle [position, mass, material] to orkan domain
--
function Orkan.addParticle(domain, position, mass, material)
	local n = #domain.x + 1
	domain.F[n] = Math.IdentityMatrix(domain.dim)
	domain.ms[n] = mass
	domain.mt[n] = material
	domain.x[n] = position
end

--
-- Bulk adds particles to orkan domain
--
function Orkan.bulkAddParticles(domain, positions, masses, materials)
	local x, ms, mt, F = domain.x, domain.ms, domain.mt, domain.F
	local identityMat = Math.IdentityMatrix(domain.dim)
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