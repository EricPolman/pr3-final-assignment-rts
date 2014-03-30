Possible optimizations
======================
High-Level
----------
* Smoke::Tick()
	Nuthin'?
	
* Bullet::Tick()
	Nuthin'?
	
* Tank::Fire()
	Keep track of next bullet position. O(n*k) to O(n)
	
* Tank::Tick()
	Nuthin?
	
* Game::DrawTanks()
	Optimize glow: O(n*k) where K=1024*768, could possibly be K=40*40
	
* Game::Tick()
	

Low-Level
---------
* Smoke::Tick()
	Get rid of !(frame++ & 7), put in constructor of smoke
	
* Bullet::Tick()
	
* Tank::Fire()
	Send tank's flags, saves an OR
	
* Tank::Tick()
	
* Game::DrawTanks()
	Optimize glow: Precalculate distances
	Float -> Int conversions, ew.
	
* Game::Tick()
	Drawing backdrop: optimize CopyTo (one memset.)
	Update bullets: Early exit
	
SIMD
----
* Smoke::Tick()
	4 puffs at the same time?
	
* Bullet::Tick()
	4 bullets at the same time? Lots of masking
	
* Tank::Fire()
	
	
* Tank::Tick()
	4 mountains at the same time?
	
* Game::DrawTanks()
	
* Game::Tick()
	

More
----
Tanks ** to Tanks *
Save opponent start variable
Fix mouse lag