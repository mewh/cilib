/*
 * Copyright (C) 2003 - 2009
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
package net.sourceforge.cilib.io.transform;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import net.sourceforge.cilib.io.DataTable;
import net.sourceforge.cilib.io.StandardDataTable;
import net.sourceforge.cilib.io.exception.CIlibIOException;
import net.sourceforge.cilib.type.types.Bit;
import net.sourceforge.cilib.type.types.Real;
import net.sourceforge.cilib.type.types.StringType;
import net.sourceforge.cilib.type.types.Type;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 * Test the Type conversion operator. Note the untyped usage of DataTable, the
 * type is infered.
 * @author andrich
 */
public class TypeConversionOperatorTest {

    private DataTable dataTable;
    private ArrayList<List<Type>> referenceRows;

    @Before
    public void setupTestData() {
        referenceRows = new ArrayList<List<Type>>();
        dataTable = new StandardDataTable<StringType>();

        List<String> tokens = Arrays.asList("0","0.2","0.3E-6","TRUE","-0.5E-2","class1");
        dataTable.addRow(tokens);
        List<Type> row = new ArrayList<Type>();
        row.add(new Real(new Double(tokens.get(0))));
        row.add(new Real(new Double(tokens.get(1))));
        row.add(new Real(new Double(tokens.get(2))));
        row.add(new Bit(true));
        row.add(new Real(new Double(tokens.get(4))));
        row.add(new StringType(tokens.get(5)));
        referenceRows.add(row);

        tokens = Arrays.asList("1","1.2","1.3E-6","T","-1.5E-3","class2");
        dataTable.addRow(tokens);
        row = new ArrayList<Type>();
        row.add(new Real(new Double(tokens.get(0))));
        row.add(new Real(new Double(tokens.get(1))));
        row.add(new Real(new Double(tokens.get(2))));
        row.add(new Bit(true));
        row.add(new Real(new Double(tokens.get(4))));
        row.add(new StringType(tokens.get(5)));
        referenceRows.add(row);

        tokens = Arrays.asList("2","2.2","2.3E-2","false","-2.5E-2","class3");
        dataTable.addRow(tokens);
        row = new ArrayList<Type>();
        row.add(new Real(new Double(tokens.get(0))));
        row.add(new Real(new Double(tokens.get(1))));
        row.add(new Real(new Double(tokens.get(2))));
        row.add(new Bit(false));
        row.add(new Real(new Double(tokens.get(4))));
        row.add(new StringType(tokens.get(5)));
        referenceRows.add(row);
    }

    @Test
    public void testOperator() throws CIlibIOException {
        TypeConversionOperator operator = new TypeConversionOperator();
        DataTable resultTable = operator.operate(dataTable);

        int size = referenceRows.size();

        for (int i = 0; i < size; i++) {
            Assert.assertEquals(referenceRows.get(i), resultTable.getRow(i));
        }
    }
}
