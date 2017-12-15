module F2 where

import Data.List

data SeqType = DNA|Protein deriving (Show, Eq)
data MolSeq = MolSeq { seqName :: String, seqSequence :: String, seqType :: SeqType} deriving (Show) 
type Matrix = [[(Char, Int)]] 
data Profile = Profile {profileName :: String, profileMatrix :: Matrix, numSeq :: Int, profileType :: SeqType} deriving (Show)

class Evol a where
	distance :: a -> a -> Double 
	name :: a -> String 
	distanceMatrix :: [a] -> [(String,String,Double)] 
	distanceMatrix [] = [] 
	distanceMatrix list = (calcDist (head list) list) ++ (distanceMatrix (tail list))
		where 
			calcDist x [] = [] 
			calcDist x (y:xs) = (name x, name y, distance x y):calcDist x xs

instance Evol MolSeq where
	distance a b = seqDistance a b 
	name a = seqName a 

instance Evol Profile where
	distance a b = profileDistance a b 
	name a = profileName a

nucleotides = "ACGT"
aminoacids = sort "ARNDCEQGHILKMFPSTWYVX"

isDNA :: String -> Bool 
isDNA [] = True 
isDNA (x:xs)
	| x `elem` nucleotides = isDNA xs
	| otherwise = False

string2seq :: String -> String -> MolSeq
string2seq name sq 
	| isDNA sq = MolSeq name sq DNA
	| otherwise = MolSeq name sq Protein

seqLength :: MolSeq -> Int
seqLength mol = length (seqSequence mol) 

isSameType :: MolSeq -> MolSeq -> Bool
isSameType a b = seqType a == seqType b

seqDistance :: MolSeq -> MolSeq -> Double
seqDistance seqA seqB 
	| not (isSameType seqA seqB)  = error "Cannot compare DNA with Protein sequence."
	| seqType seqA == DNA = dnaDist (seqSequence seqA) (seqSequence seqB)
	| otherwise = proteinDist (seqSequence seqA) (seqSequence seqB)
	where
		dnaDist s1 s2 
			| a > 0.74 = 3.3  
			| otherwise = -(3/4)*log(1-(4*a)/3)
		proteinDist s1 s2
			| a >= 0.94 = 3.7
			| otherwise = -(19/20)*log(1-(20*a)/19)
		a = normHammingDist (seqSequence seqA) (seqSequence seqB) 0 0 
			where 
				normHammingDist [] sq2 diff sqLength = let extra = length sq2 in (fromIntegral(diff+extra))/(fromIntegral(sqLength+extra))
				normHammingDist sq1 [] diff sqLength = let extra = length sq1 in (fromIntegral(diff+extra))/(fromIntegral(sqLength+extra))
				normHammingDist (s1:sq1) (s2:sq2) diff sqLength 
					| s1 == s2 = normHammingDist sq1 sq2 diff (sqLength+1) 
					| otherwise = normHammingDist sq1 sq2 (diff+1) (sqLength+1)


molseqs2profile :: String -> [MolSeq] -> Profile
molseqs2profile name seqList = Profile {profileName=name, profileMatrix=makeProfileMatrix seqList, numSeq=length seqList, profileType=seqType (head seqList)}

makeProfileMatrix :: [MolSeq] -> Matrix
makeProfileMatrix [] = error "Empty sequence list" 
makeProfileMatrix sl = res
  where 
    t = seqType (head sl) 
    defaults = 
      if (t == DNA) then
      	
        zip nucleotides (replicate (length nucleotides) 0) 
      else 
        zip aminoacids (replicate (length aminoacids) 0)   
    strs = map seqSequence sl                              
    tmp1 = map (map (\x -> ((head x), (length x))) . group . sort) (transpose strs)
    equalFst a b = (fst a) == (fst b)
    res = map sort (map (\l -> unionBy equalFst l defaults) tmp1)


profileFrequency :: Profile -> Int -> Char -> Double

profileFrequency p i c = (fromIntegral(findCount (matrix !! i)))/(fromIntegral(numSeq p)) 
	where 
		matrix = profileMatrix p  
		findCount [] = 0 
		findCount (h:t) = if (fst h) == c 
							
							then (snd h)  
							else findCount t  

profileDistance :: Profile -> Profile -> Double
profileDistance p1 p2 
	| profileType p1 == DNA = loop (length (profileMatrix p1)) True 0 
	| otherwise = loop (length (profileMatrix p1)) False 0
	where
		loop n dna res 
			| n<=0 = res 
			| dna = col nucleotides (n-1) dna res 
			| otherwise = col aminoacids (n-1) dna res

		col [] n dna res = loop n dna res
		col (c:chars) n dna res = col chars n dna (res+abs((profileFrequency p1 n c)-(profileFrequency p2 n c)))







