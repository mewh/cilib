/**
 * Copyright (C) 2003 - 2008
 * Computational Intelligence Research Group (CIRG@UP)
 * Department of Computer Science
 * University of Pretoria
 * South Africa
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.sourceforge.cilib.pso.niching;

import java.util.HashMap;
import java.util.Map;
import net.sourceforge.cilib.container.Pair;
import net.sourceforge.cilib.entity.Entity;
import net.sourceforge.cilib.problem.Fitness;

/**
 *
 * @author gpampara
 */
public class StandardSwarmCreationStrategy implements SwarmCreationStrategy {

	private Map<Entity, Pair<Fitness, Integer>> deviations;

	public StandardSwarmCreationStrategy() {
		this.deviations = new HashMap<Entity, Pair<Fitness, Integer>>();
	}

	@Override
	public void create(NichePSO algorithm) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

}