use ark_ec::models::{short_weierstrass_jacobian::GroupAffine as SWAffine, SWModelParameters};
use ark_ec::PairingEngine;
use ark_ff::{Field, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly_commit::marlin::marlin_pc::MarlinKZG10;
use ark_poly_commit::{Polynomial, PolynomialCommitment};



pub(crate) enum Node<P, C, V, const DEPTH: usize, const WIDTH: usize>
where
    P: SWModelParameters,
    C: CommitmentScheme<P>,
    V: ToField<P> + Debug,
    Fr<P>: From<i64>,
{
    Internal {
        commitment: C::Commitment,
        children: HashMap<[u8; 4], Self>,
    },
    Value {
        commitment: C::Commitment,
        stem: BitVec,
        values_commitments: [C::Commitment; 1],
        values: HashMap<[u8; 4], V>,
    },
}




impl<F: Field, P: Polynomial<F>, PC: PolynomialCommitment<F, P>> VerkleTree<F, P, PC>
where
    PC::Commitment: ToFieldElements<F>,
{
    /// Create a verkle tree with the given depth and branching factor
    pub fn new(comm_key: PC::CommitterKey, depth: usize, branching_factor: usize) -> Self {
        ///New
        let num_limbs = L::num_limbs();
        assert!(
            committer.get_domain_size() >= num_limbs + 2,
            "domain size must be larger than or equal to {}",
            num_limbs + 2
        );
        Self {
            root: VerkleNode::default(),
            committer,
        }
    }

    pub fn get_width(&self) -> usize {
        self.committer.get_domain_size()
    }

    /// Returns the depth of the tree
    pub fn depth(&self) -> usize {
        self.root.path::len
    }

    /// Returns the polynomial commitment at the root of the tree
    pub fn root(&self) -> PC::Commitment {
        self.root::PC::Commitment
    }

    /// Add an element to the tree at the given position
    pub fn insert(&mut self, position: usize, x: F) {
        self.root.insert(key.to_path(), key, value)
    }

    /// Batch-open the verkle tree at the given set of positions
    pub fn open(&self, position: Vec<usize>) -> Option<(Vec<F>, VerkleProof<F, P, PC>)> {
        self.root.get(key.to_path(), key)
    }

    /// Check the correctness of an opening
    pub fn check(
        root: PC::Commitment,
        vk: PC::VerifierKey,
        (x, proof): (Vec<F>, VerkleProof<F, P, PC>),
    ) -> bool {
        
    }
}

pub trait ToFieldElements<F: Field> {
    //  a possible method for converting a polynomial commitment into an vector of field
    // elements.
    fn to_field_elements(&self) -> Vec<F>;
}

impl<'a, 'b, P, E> ToFieldElements<P::ScalarField>
    for ark_poly_commit::marlin::marlin_pc::Commitment<E>
where
    P: SWModelParameters,
    E: PairingEngine<Fq = P::BaseField, G1Affine = SWAffine<P>>,
    P::ScalarField: PrimeField,
    P::BaseField: PrimeField<BigInt = <P::ScalarField as PrimeField>::BigInt>,
{
    fn to_field_elements(&self) -> Vec<P::ScalarField> {
        // unused degree bounds, and so ignore the shifted part
        let _ = self.shifted_comm;
        [self.comm.0.x, self.comm.0.y]
            .iter()
            .map(|a| P::ScalarField::from_repr(a.into_repr()).unwrap())
            .collect()
    }
}


enum VerkleTree<F: Field, P: Polynomial<F>, PC: PolynomialCommitment<F, P>> {
    K: AbstractKey,
    L: LeafNodeValue<K, GA>,
    GA: CurveAffine,
{
    Leaf {
        path: K::Path,
        stem: K::Stem,
        info: L,
    },
    Internal {
        path: K::Path,
        children: HashMap<usize, VerkleNode<K, L, GA>>, // HashMap<u8, VerkleNode<K, V, GA>>
        info: InternalNodeValue<GA>,
    },
}

}

struct VerkleProof<F: Field, P: Polynomial<F>, PC: PolynomialCommitment<F, P>> {
    // Just here to get rid of the unused variable warning
    P: SWModelParameters,
    C: CommitmentScheme<P>,
    V: ToField<P> + Clone + Debug,
    Fr<P>: From<i64>,
{
    key: BitVec,
    value: V,
    path: Vec<(C::Commitment, C::Opening)>,
    value_commitment: C::Commitment,
    value_opening: C::Opening,
}

use ark_bn254::{Bn254, Fr};
type Bn254KZG = MarlinKZG10<Bn254, DensePolynomial<Fr>>;

pub trait AbstractKey: Clone + Copy + Debug + Ord {
    type Stem: AbstractStem + Default;
    type Path: AbstractPath + Default;

    fn get_stem(&self) -> Self::Stem;

    fn get_suffix(&self) -> usize;

    fn to_path(&self) -> Self::Path;

    fn get_branch_at(&self, depth: usize) -> usize;
}

pub trait AbstractValue: Clone + Copy + Debug + Eq {}

pub trait AbstractPath: Sized + Clone + Debug + Default + Eq + IntoIterator<Item = usize> {
    type RemovePrefixError: Debug + Send + Sync + 'static;

    fn get_next_path_and_branch(&self) -> (Self, usize);

    fn get_suffix(&self) -> usize;

    fn is_proper_prefix_of(&self, full_path: &Self) -> bool;

    fn remove_prefix(&self, prefix: &Self) -> Result<Self, Self::RemovePrefixError>;

    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn push(&mut self, value: usize);
}

pub trait AbstractStem: Clone + Debug + Eq {
    type Path: AbstractPath;

    fn to_path(&self) -> Self::Path;
}

pub trait IntoFieldElement<F: PrimeField> {
    type Err: Send + Sync + 'static;

    fn into_field_element(self) -> Result<F, Self::Err>;
}

pub trait NodeValue<GA>
where
    GA: CurveAffine,
{
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn get_digest_mut(&mut self) -> &mut Option<GA::Scalar>;

    fn get_digest(&self) -> Option<&GA::Scalar>;
}

pub trait LeafNodeValue<K, GA>: Clone + Default + NodeValue<GA>
where
    K: AbstractKey,
    GA: CurveAffine,
{
    type Value: AbstractValue;

    fn new() -> Self;

    fn insert(&mut self, key: usize, value: Self::Value) -> Option<Self::Value>;

    fn get(&self, key: &usize) -> Option<&Self::Value>;

    fn remove(&mut self, key: &usize) -> Option<Self::Value>;

    fn compute_digest<C: Committer<GA>>(
        &mut self,
        stem: &mut K::Stem,
        committer: &C,
    ) -> anyhow::Result<GA::Scalar>;

    
    fn bits_of_value() -> usize;

    fn num_limbs() -> usize;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InternalNodeValue<GA>
where
    GA: CurveAffine,
{
    /// The number of children which are `Some` rather than `None`.
    num_nonempty_children: usize,

    /// The commitment of this node.
    /// If it has not computed yet, `commitment` set `None`.
    commitment: Option<GA>,

    /// The digest of `commitment`.
    /// If it has not computed yet, `digest` set `None`.
    digest: Option<GA::Scalar>,
}

impl<GA> NodeValue<GA> for InternalNodeValue<GA>
where
    GA: CurveAffine,
{
    fn len(&self) -> usize {
        self.num_nonempty_children
    }

    fn get_digest_mut(&mut self) -> &mut Option<GA::Scalar> {
        &mut self.digest
    }

    fn get_digest(&self) -> Option<&GA::Scalar> {
        let digest = &self.digest;

        digest.into()
    }
}

impl<GA> InternalNodeValue<GA>
where
    GA: CurveAffine,
{
    pub fn get_commitment_mut(&mut self) -> &mut Option<GA> {
        &mut self.commitment
    }

    pub fn get_commitment(&self) -> Option<&GA> {
        let commitment = &self.commitment;

        commitment.into()
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum VerkleNode<K, L, GA>
where
    K: AbstractKey,
    L: LeafNodeValue<K, GA>,
    GA: CurveAffine,
{
    Leaf {
        path: K::Path,
        stem: K::Stem,
        info: L,
    },
    Internal {
        path: K::Path,
        children: HashMap<usize, VerkleNode<K, L, GA>>, // HashMap<u8, VerkleNode<K, V, GA>>
        info: InternalNodeValue<GA>,
    },
}

impl<K, L, GA> VerkleNode<K, L, GA>
where
    K: AbstractKey,
    L: LeafNodeValue<K, GA>,
    GA: CurveAffine,
    GA::Base: PrimeField,
{
    pub fn new_leaf_node_with_entry(path: K::Path, key: K, value: L::Value) -> Self {
        // let mut leaves = HashMap::new();
        // leaves.insert(key.get_suffix(), value);
        let mut info = L::new();
        info.insert(key.get_suffix(), value);
        Self::Leaf {
            stem: key.get_stem(),
            path,
            info,
        }
    }

    pub fn new_internal_node_with_children(
        path: K::Path,
        children: HashMap<usize, VerkleNode<K, L, GA>>,
    ) -> Self {
        let num_nonempty_children = children.len();
        Self::Internal {
            path,
            children,
            info: InternalNodeValue {
                num_nonempty_children,
                commitment: None,
                digest: None,
            },
        }
    }

    pub fn get_path(&self) -> &K::Path {
        match self {
            Self::Leaf { path, .. } => path,
            Self::Internal { path, .. } => path,
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            Self::Leaf { info, .. } => info.is_empty(),
            Self::Internal { info, .. } => info.is_empty(),
        }
    }

    

    pub fn get_digest(&self) -> Option<&GA::Scalar> {
        match self {
            Self::Leaf { info, .. } => info.get_digest(),
            Self::Internal { info, .. } => info.get_digest(),
        }
    }
}

impl<K, L, GA> Default for VerkleNode<K, L, GA>
where
    K: AbstractKey,
    L: LeafNodeValue<K, GA>,
    GA: CurveAffine,
    GA::Base: PrimeField,
{
    fn default() -> Self {
        Self::new_internal_node_with_children(K::Path::default(), HashMap::new())
    }
}

impl<P, K, L, GA> VerkleNode<K, L, GA>
where
    P: Default + AbstractPath,
    K: AbstractKey<Path = P>,
    K::Stem: AbstractStem<Path = P> + IntoFieldElement<GA::Scalar>,
    L: LeafNodeValue<K, GA>,
    GA: CurveAffine,
    <GA as CurveAffine>::Base: PrimeField,
{
    pub fn insert(&mut self, relative_path: K::Path, key: K, value: L::Value) -> Option<L::Value> {
        if relative_path.is_empty() {
            panic!("`relative_path` must be non-empty.");
        }

        match self {
            VerkleNode::Leaf {
                stem, path, info, ..
            } => {
                let stem_relative_path = stem
                    .to_path()
                    .remove_prefix(path)
                    .expect("unreachable code");
                if stem_relative_path.is_empty() {
                    panic!("`relative_path` must be non-empty.");
                }
                if key.get_stem().eq(stem) {
                    return info.insert(key.get_suffix(), value);
                }

                // A new branch node has to be inserted. Depending
                // on the next branch in both keys, a recursion into
                // the moved leaf node can occur.
                let (_, next_branch_of_existing_key) =
                    stem_relative_path.get_next_path_and_branch();
                // assert!(next_branch_of_existing_key < WIDTH);

                let mut new_branch = {
                    let mut children = HashMap::new();
                    let mut new_path = path.clone();
                    new_path.push(next_branch_of_existing_key);
                    let moving_child = VerkleNode::Leaf {
                        stem: stem.clone(),
                        path: new_path,
                        info: info.clone(),
                    };
                    children.insert(next_branch_of_existing_key, moving_child);

                    VerkleNode::new_internal_node_with_children(path.clone(), children)
                };

                let (_, next_branch_of_inserting_key) = relative_path.get_next_path_and_branch();
                // assert!(next_branch_of_inserting_key < WIDTH);

                if next_branch_of_inserting_key != next_branch_of_existing_key {
                    // Next branch differs, so this was the last level.
                    // Insert it directly into its suffix.
                    let mut info = L::new();
                    info.insert(key.get_suffix(), value);
                    let mut new_path = path.clone();
                    new_path.push(next_branch_of_inserting_key);
                    let leaf_node = VerkleNode::Leaf {
                        stem: key.get_stem(),
                        path: new_path,
                        info,
                    };

                    match &mut new_branch {
                        VerkleNode::Internal { children, info, .. } => {
                            children.insert(next_branch_of_inserting_key, leaf_node);
                            info.num_nonempty_children += 1;
                        }
                        VerkleNode::Leaf { .. } => {
                            panic!("unreachable code");
                        }
                    }
                    let _ = std::mem::replace(self, new_branch);

                    return None;
                }

                let _ = std::mem::replace(self, new_branch);

                self.insert(relative_path, key, value)
            }
            VerkleNode::Internal {
                path,
                children,
                info,
                ..
            } => {
                let _ = info.commitment.take();
                let _ = info.digest.take();

                let (next_relative_path, next_branch_of_inserting_key) =
                    relative_path.get_next_path_and_branch();
                // assert!(next_branch_of_inserting_key < WIDTH);

                if let Some(child) = children.get_mut(&next_branch_of_inserting_key) {
                    child.insert(next_relative_path, key, value)
                } else {
                    let mut new_path = path.clone();
                    new_path.push(next_branch_of_inserting_key);
                    children.insert(
                        next_branch_of_inserting_key,
                        VerkleNode::new_leaf_node_with_entry(new_path, key, value),
                    );

                    None
                }
            }
        }
    }

    pub fn remove(&mut self, relative_path: K::Path, key: &K) -> Option<L::Value> {
        match self {
            Self::Leaf { info, .. } => info.remove(&key.borrow().get_suffix()),
            Self::Internal {
                children,
                info:
                    InternalNodeValue {
                        commitment, digest, ..
                    },
                ..
            } => {
                let _ = commitment.take();
                let _ = digest.take();

                let (next_path, next_branch) = relative_path.get_next_path_and_branch();
                // assert!(next_branch < WIDTH);

                if let Some(child) = children.get_mut(&next_branch) {
                    let old_value = child.remove(next_path, key);

                    // Remove a empty node if any.
                    if child.is_empty() {
                        let _ = children.remove(&next_branch);
                    }

                    old_value
                } else {
                    None
                }
            }
        }
    }

    /// Get a value from this tree.
    pub fn get(&self, relative_path: K::Path, key: &K) -> Option<&L::Value> {
        match &self {
            Self::Leaf { stem, info, .. } => {
                if key.get_stem() != stem.clone() {
                    None
                } else {
                    info.get(&key.get_suffix())
                }
            }
            Self::Internal { children, .. } => {
                let (next_path, next_branch) = relative_path.get_next_path_and_branch();
                // assert!(next_branch < WIDTH);

                children
                    .get(&next_branch)
                    .and_then(|child| child.get(next_path, key))
            }
        }
    }


}

pub fn compute_commitment_of_internal_node<GA: CurveAffine, C: Committer<GA>>(
    committer: &C,
    children_digests: Vec<GA::Scalar>,
) -> anyhow::Result<GA> {
    committer
        .commit(&children_digests)
        .or_else(|_| anyhow::bail!("Fail to make a commitment of given polynomial."))
}

impl<P, K, L, GA> VerkleNode<K, L, GA>
where
    P: Default + AbstractPath,
    K: AbstractKey<Path = P>,
    K::Stem: AbstractStem<Path = P> + IntoFieldElement<GA::Scalar>,
    L: LeafNodeValue<K, GA>,
    GA: CurveAffine,
    GA::Base: PrimeField,
{
    pub fn compute_digest<C: Committer<GA>>(
        &mut self,
        committer: &C,
    ) -> anyhow::Result<GA::Scalar> {
        if let Some(d) = self.get_digest() {
            return Ok(*d);
        }

        match self {
            VerkleNode::Leaf { stem, info, .. } => info.compute_digest::<C>(stem, committer),
            VerkleNode::Internal { children, info, .. } => {
                let width = committer.get_domain_size();
                let mut children_digests = vec![GA::Scalar::zero(); width];
                for (&i, child) in children.iter_mut() {
                    children_digests[i] = child.compute_digest(committer)?;
                }

                let tmp_commitment =
                    compute_commitment_of_internal_node(committer, children_digests)?;
                let tmp_digest = point_to_field_element(&tmp_commitment)?;

                let _ = std::mem::replace(&mut info.commitment, Some(tmp_commitment));
                let _ = std::mem::replace(&mut info.digest, Some(tmp_digest));

                Ok(tmp_digest)
            }
        }
    }
}
