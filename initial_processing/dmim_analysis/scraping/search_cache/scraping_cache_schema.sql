-- table for caching PubChem data
CREATE TABLE pubchem (
    -- PubChem compound ID (CID)
    cid INTEGER UNIQUE NOT NULL,
    -- compound name
    name TEXT NOT NULL,
    -- monoisotopic mass
    monoiso REAL,
    -- SMILES structure
    smiles TEXT
);

-- table for caching images from
CREATE TABLE cactus (
    -- unique integer identifier
    i INTEGER UNIQUE NOT NULL,
    -- SMILES (or InChI) structure
    structure TEXT NOT NULL,
    -- png image of the structure
    img BLOB NOT NULL
);
